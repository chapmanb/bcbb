"""SeqIO compatible functionality for Generic Feature Format (GFF) files.

This deals with GFF3 formatted files, a tab delimited format for storing
sequence features and annotations:

http://www.sequenceontology.org/gff3.shtml
"""
import collections
from Bio.SeqFeature import SeqFeature, FeatureLocation

class GFFFeatureAdder:
    """Move through a GFF file, adding new features to SeqRecord objects.

    This class is instantiated with a dictionary where the keys are IDs
    corresponding to those in the first column of the GFF file. As the GFF file
    is processed, the items are added to the appropriate record as features.
    """
    def __init__(self, base_dict):
        self.base = base_dict
        self._strand_map = {'+' : 1, '-' : -1, '?' : None, None: None}
        self._filter_info = dict(gff_id = [0], gff_types = [1, 2])

    def limit_names(self):
        return self._filter_info.keys()

    def available_limits(self, gff_handle):
        """Return dictionary information on possible limits for this file.

        This returns a nested dictionary with the following structure:
        
        keys -- names of items to filter by
        values -- dictionary with:
            keys -- filter choice
            value -- counts of that filter in this file
        """
        cur_limits = dict()
        for filter_key in self._filter_info.keys():
            cur_limits[filter_key] = collections.defaultdict(int)
        for line in gff_handle:
            # ignore comment lines
            if line.strip()[0] != "#":
                parts = [p.strip() for p in line.split('\t')]
                assert len(parts) == 9, line
                for filter_key, cur_indexes in self._filter_info.items():
                    cur_id = tuple([parts[i] for i in cur_indexes])
                    cur_limits[filter_key][cur_id] += 1
        # get rid of the default dicts
        final_dict = dict()
        for key, value_dict in cur_limits.items():
            if len(key) == 1:
                key = key[0]
            final_dict[key] = dict(value_dict)
        return final_dict

    def add_features(self, gff_handle, limit_info=None):
        """Add features to the base records from the current GFF handle.

        limit_info is a dictionary of items to help limit which items are
        added:
            keys -- limit_name specifying an item to limit by
            value -- list of items to filter by
        The available_limits function returns a dictionary of limits and
        counts for querying new GFF files.
        """
        parents = collections.defaultdict(list)
        children = collections.defaultdict(list)
        for line in gff_handle:
            # ignore comment lines XXX We could handle directives (##) here
            if line.strip()[0] != "#":
                parts = [p.strip() for p in line.split('\t')]
                assert len(parts) == 9, line
                should_do = True
                if limit_info:
                    for limit_name, limit_values in limit_info.items():
                        cur_id = tuple([parts[i] for i in
                            self._filter_info[limit_name]])
                        if len(cur_id) == 1:
                            cur_id = cur_id[0]
                        if cur_id not in limit_values:
                            should_do = False
                            break
                if should_do:
                    rec = self.base[parts[0]]
                    try:
                        rec, parents, children = self._add_gff_line(rec, parts,
                                parents, children)
                    except KeyboardInterrupt:
                        print line
                        raise
        self._add_parent_child_features(dict(parents), dict(children))

    def _add_parent_child_features(self, parents, children):
        """Add nested features with parent child relationships.
        """
        for rec_id, all_parents in parents.items():
            rec = self.base[rec_id]
            for cur_parent in all_parents:
                cur_parent, children = self._add_children_to_parent(cur_parent,
                        children)
                rec.features.append(cur_parent)
        if len(children) > 0:
            raise ValueError("Non-added children: %s" % children.keys())

    def _add_children_to_parent(self, cur_parent, children):
        """Recursively add children to parent features.
        """
        if children.has_key(cur_parent.id):
            cur_children = children[cur_parent.id]
            for cur_child in cur_children:
                cur_child, children = self._add_children_to_parent(cur_child,
                        children)
                cur_parent.sub_features.append(cur_child)
            del children[cur_parent.id]
        return cur_parent, children

    def _add_gff_line(self, rec, gff_parts, parents, children):
        """Add details from a GFF line to the given SeqRecord.
        """
        gff_parts = [(None if p == '.' else p) for p in gff_parts]
        assert rec.id == gff_parts[0], "ID mismatch: %s %s" % (rec.id,
                gff_parts[0])
        # collect all of the base qualifiers for this item
        quals = collections.defaultdict(list)
        if gff_parts[1]:
            quals["source"].append(gff_parts[1])
        if gff_parts[5]:
            quals["score"].append(gff_parts[5])
        if gff_parts[7]:
            quals["phase"].append(gff_parts[7])
        for key, val in [a.split('=') for a in gff_parts[8].split(';')]:
            quals[key].extend(val.split(','))
        quals = dict(quals)
        # if we are describing a location, then we are a feature
        if gff_parts[3] and gff_parts[4]:
            #if quals.has_key('ID') or quals.has_key('Parent'):
            #    print gff_parts[1:6], quals
            location = FeatureLocation(int(gff_parts[3]) - 1, int(gff_parts[4]))
            new_feature = SeqFeature(location, gff_parts[2],
                    id = quals.get('ID', [''])[0],
                    strand = self._strand_map[gff_parts[6]])
            new_feature.qualifiers = quals
            # Handle flat features
            if not new_feature.id:
                rec.features.append(new_feature)
            # features that have parents need to link so we can pick up
            # the relationship
            elif new_feature.qualifiers.has_key('Parent'):
                for parent in new_feature.qualifiers['Parent']:
                    children[parent].append(new_feature)
            # top level features
            else:
                parents[rec.id].append(new_feature)
        # otherwise, associate these annotations with the full record
        else:
            # add these as a list of annotations, checking not to overwrite
            # current values
            for key, vals in quals:
                if rec.annotations.has_key(key):
                    try:
                        rec.annotations[key].extend(vals)
                    except AttributeError:
                        rec.annotations[key] = [rec.annotations[key]] + vals
                else:
                    rec.annotations[key] = vals
        return rec, parents, children
