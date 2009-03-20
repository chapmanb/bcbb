"""Parse GFF files into features attached to Biopython SeqRecord objects.

This deals with GFF3 formatted files, a tab delimited format for storing
sequence features and annotations:

http://www.sequenceontology.org/gff3.shtml

The implementation utilizes map/reduce parsing of GFF using Disco. Disco
(http://discoproject.org) is a Map-Reduce framework for Python utilizing
Erlang for parallelization. The code works on a single processor without
Disco using the same architecture.
"""
import os
import collections

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def gff_line_map(line, params):
    """Map part of Map-Reduce; parses a line of GFF into a dictionary.
    """
    strand_map = {'+' : 1, '-' : -1, '?' : None, None: None}
    line = line.strip()
    if line[0] != "#":
        parts = line.split('\t')
        should_do = True
        if params.limit_info:
            for limit_name, limit_values in params.limit_info.items():
                cur_id = tuple([parts[i] for i in 
                    params.filter_info[limit_name]])
                if cur_id not in limit_values:
                    should_do = False
                    break
        if should_do:
            assert len(parts) == 9, line
            gff_parts = [(None if p == '.' else p) for p in parts]
            gff_info = dict()
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
            gff_info['quals'] = dict(quals)
            # if we are describing a location, then we are a feature
            if gff_parts[3] and gff_parts[4]:
                gff_info['location'] = [int(gff_parts[3]) - 1,
                        int(gff_parts[4])]
                gff_info['type'] = gff_parts[2]
                gff_info['id'] = quals.get('ID', [''])[0]
                gff_info['rec_id'] = gff_parts[0]
                gff_info['strand'] = strand_map[gff_parts[6]]
                # Handle flat features
                if not gff_info['id']:
                    final_key = 'feature'
                # features that have parents need to link so we can pick up
                # the relationship
                elif gff_info['quals'].has_key('Parent'):
                    final_key = 'child'
                # top level features
                else:
                    final_key = 'parent'
            # otherwise, associate these annotations with the full record
            else:
                final_key = 'annotation'
            return [(final_key, (simplejson.dumps(gff_info) if params.jsonify
                else gff_info))]
    return []

def gff_line_reduce(map_results, out, params):
    """Reduce part of Map-Reduce; combines results of parsed features.
    """
    final_items = dict()
    for gff_type, final_val in map_results:
        send_val = (simplejson.loads(final_val) if params.jsonify else 
                final_val)
        try:
            final_items[gff_type].append(send_val)
        except KeyError:
            final_items[gff_type] = [send_val]
    for key, vals in final_items.items():
        out.add(key, (simplejson.dumps(vals) if params.jsonify else vals))

class GFFMapReduceFeatureAdder:
    """Move through a GFF file, adding new features to SeqRecord objects.
    """
    def __init__(self, base_dict, disco_host=None, create_missing=True):
        """Initialize with dictionary of records to add to.

        This class is instantiated with a dictionary where the keys are IDs
        corresponding to those in the first column of the GFF file. As the GFF
        file is processed, the items are added to the appropriate record as
        features.

        disco_host - Web reference to a Disco (http://discoproject.org) host
        which will be used for parallelizing the GFF reading job.

        create_missing - If True, create blank records for GFF ids not in
        the base_dict. If False, an error will be raised.
        """
        self.base = base_dict
        self._create_missing = create_missing
        self._map_fn = gff_line_map
        self._reduce_fn = gff_line_reduce
        self._disco_host = disco_host
        # details on what we can filter items with
        self._filter_info = dict(gff_id = [0], gff_types = [1, 2])

    def add_features(self, gff_file, limit_info=None):
        # turn all limit information into tuples for identical comparisons
        final_limit_info = {}
        for key, values in limit_info.items():
            final_limit_info[key] = [tuple(v) for v in values]
        if self._disco_host:
            results = self._disco_process(gff_file, final_limit_info)
        else:
            results = self._std_process(gff_file, final_limit_info)
        #print results
        self._add_annotations(results.get('annotation', []))
        [self._add_toplevel_feature(f) for f in results.get('feature', [])]
        self._add_parent_child_features(results.get('parent', []),
                results.get('child', []))
    
    def _add_parent_child_features(self, parents, children):
        """Add nested features with parent child relationships.
        """
        children_prep = collections.defaultdict(list)
        for child_dict in children:
            child_feature = self._get_feature(child_dict)
            for parent in child_feature.qualifiers['Parent']:
                children_prep[parent].append(child_feature)
        children = dict(children_prep)
        for cur_parent_dict in parents:
            cur_parent = self._add_toplevel_feature(cur_parent_dict)
            cur_parent, children = self._add_children_to_parent(cur_parent,
                    children)
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

    def _add_annotations(self, anns):
        """Add annotation data from the GFF file to records.
        """
        # add these as a list of annotations, checking not to overwrite
        # current values
        for ann in anns:
            rec = self.base[ann['rec_id']]
            for key, vals in ann['quals']:
                if rec.annotations.has_key(key):
                    try:
                        rec.annotations[key].extend(vals)
                    except AttributeError:
                        rec.annotations[key] = [rec.annotations[key]] + vals
                else:
                    rec.annotations[key] = vals

    def _add_toplevel_feature(self, feature_dict):
        """Add a toplevel non-nested feature to the appropriate record.
        """
        new_feature = self._get_feature(feature_dict)
        try:
            self.base[feature_dict['rec_id']].features.append(new_feature)
        except KeyError:
            if not self._create_missing:
                raise
            else:
                new_rec = SeqRecord(Seq(""), feature_dict['rec_id'])
                new_rec.features.append(new_feature)
                self.base[feature_dict['rec_id']] = new_rec
        return new_feature

    def _get_feature(self, feature_dict):
        """Retrieve a Biopython feature from our dictionary representation.
        """
        location = FeatureLocation(*feature_dict['location'])
        new_feature = SeqFeature(location, feature_dict['type'],
                id=feature_dict['id'], strand=feature_dict['strand'])
        new_feature.qualifiers = feature_dict['quals']
        return new_feature

    def _std_process(self, gff_file, limit_info):
        """Process GFF addition without any parallelization.
        """
        class _LocalParams:
            def __init__(self):
                self.jsonify = False
        params = _LocalParams()
        params.limit_info = limit_info
        params.filter_info = self._filter_info
        class _LocalOut:
            def __init__(self):
                self._items = dict()
            def add(self, key, vals):
                try:
                    self._items[key].extend(vals)
                except KeyError:
                    self._items[key] = vals
            def get_results(self):
                return self._items
        out_info = _LocalOut()
        in_handle = open(gff_file)
        for line in in_handle:
            results = self._map_fn(line, params)
            self._reduce_fn(results, out_info, params)
        in_handle.close()
        return out_info.get_results()

    def _disco_process(self, gff_file, limit_info):
        """Process GFF addition, using Disco to parallelize the process.
        """
        # make these imports local; only need them when using disco
        import simplejson
        import disco
        full_file = os.path.join(os.getcwd(), gff_file)
        results = disco.job(self._disco_host, name="gff_reader",
                input=[full_file],
                params=disco.Params(limit_info=limit_info, jsonify=True,
                    filter_info=self._filter_info),
                required_modules=["simplejson", "collections"],
                map=self._map_fn, reduce=self._reduce_fn)
        processed = dict()
        for out_key, out_val in disco.result_iterator(results):
            processed[out_key] = simplejson.loads(out_val)
        return processed
