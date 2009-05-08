"""Parse GFF files into features attached to Biopython SeqRecord objects.

This deals with GFF3 formatted files, a tab delimited format for storing
sequence features and annotations:

http://www.sequenceontology.org/gff3.shtml

It will also deal with older GFF versions (GTF/GFF2):

http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
http://mblab.wustl.edu/GTF22.html

The implementation utilizes map/reduce parsing of GFF using Disco. Disco
(http://discoproject.org) is a Map-Reduce framework for Python utilizing
Erlang for parallelization. The code works on a single processor without
Disco using the same architecture.
"""
import os
import copy
import re
import collections
import urllib

from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

def _gff_line_map(line, params):
    """Map part of Map-Reduce; parses a line of GFF into a dictionary.

    Given an input line from a GFF file, this:
    - decides if the file passes our filtering limits
    - if so:
        - breaks it into component elements
        - determines the type of attribute (flat, parent, child or annotation)
        - generates a dictionary of GFF info which can be serialized as JSON
    """
    gff3_kw_pat = re.compile("\w+=")
    def _split_keyvals(keyval_str):
        """Split key-value pairs in a GFF2, GTF and GFF3 compatible way.

        GFF3 has key value pairs like:
          count=9;gene=amx-2;sequence=SAGE:aacggagccg
        GFF2 and GTF have:           
          Sequence "Y74C9A" ; Note "Clone Y74C9A; Genbank AC024206"
          name "fgenesh1_pg.C_chr_1000003"; transcriptId 869
        """
        quals = collections.defaultdict(list)
        if keyval_str is None:
            return quals
        # ensembl GTF has a stray semi-colon at the end
        if keyval_str[-1] == ';':
            keyval_str = keyval_str[:-1]
        # GFF2/GTF has a semi-colon with at least one space after it.
        # It can have spaces on both sides; wormbase does this.
        # GFF3 works with no spaces.
        # Split at the first one we can recognize as working
        parts = keyval_str.split(" ; ")
        if len(parts) == 1:
            parts = keyval_str.split("; ")
            if len(parts) == 1:
                parts = keyval_str.split(";")
        # check if we have GFF3 style key-vals (with =)
        is_gff2 = True
        if gff3_kw_pat.match(parts[0]):
            is_gff2 = False
            key_vals = [p.split('=') for p in parts]
        # otherwise, we are separated by a space with a key as the first item
        else:
            pieces = [p.strip().split(" ") for p in parts]
            key_vals = [(p[0], " ".join(p[1:])) for p in pieces]
        for key, val in key_vals:
            val = (val[1:-1] if (len(val) > 0 and val[0] == '"' 
                                 and val[-1] == '"') else val)
            if val:
                quals[key].extend(val.split(','))
            # if we don't have a value, make this a key=True/False style
            # attribute
            else:
                quals[key].append('true')
        for key, vals in quals.items():
            quals[key] = [urllib.unquote(v) for v in vals]
        return quals, is_gff2

    def _nest_gff2_features(gff_parts):
        """Provide nesting of GFF2 transcript parts with transcript IDs.

        exons and coding sequences are mapped to a parent with a transcript_id
        in GFF2. This is implemented differently at different genome centers
        and this function attempts to resolve that and map things to the GFF3
        way of doing them.
        """
        # map protein or transcript ids to a parent
        for transcript_id in ["transcript_id", "transcriptId", "proteinId"]:
            try:
                gff_parts["quals"]["Parent"] = \
                        gff_parts["quals"][transcript_id]
                break
            except KeyError:
                pass
        # case for WormBase GFF -- everything labelled as Transcript
        if gff_parts["quals"].has_key("Transcript"):
            # parent types
            if gff_parts["type"] in ["Transcript"]:
                if not gff_parts["id"]:
                    gff_parts["id"] = gff_parts["quals"]["Transcript"][0]
            # children types
            elif gff_parts["type"] in ["intron", "exon", "three_prime_UTR",
                    "coding_exon", "five_prime_UTR", "CDS", "stop_codon",
                    "start_codon"]:
                gff_parts["quals"]["Parent"] = gff_parts["quals"]["Transcript"]
        return gff_parts

    strand_map = {'+' : 1, '-' : -1, '?' : None, None: None}
    line = line.strip()
    if line[:2] == "##":
        return [('directive', line[2:])]
    elif line[0] != "#":
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
            assert len(parts) >= 9, line
            gff_parts = [(None if p == '.' else p) for p in parts]
            gff_info = dict()
            # collect all of the base qualifiers for this item
            quals, is_gff2 = _split_keyvals(gff_parts[8])
            gff_info["is_gff2"] = is_gff2
            if gff_parts[1]:
                quals["source"].append(gff_parts[1])
            if gff_parts[5]:
                quals["score"].append(gff_parts[5])
            if gff_parts[7]:
                quals["phase"].append(gff_parts[7])
            gff_info['quals'] = dict(quals)
            gff_info['rec_id'] = gff_parts[0]
            # if we are describing a location, then we are a feature
            if gff_parts[3] and gff_parts[4]:
                gff_info['location'] = [int(gff_parts[3]) - 1,
                        int(gff_parts[4])]
                gff_info['type'] = gff_parts[2]
                gff_info['id'] = quals.get('ID', [''])[0]
                gff_info['strand'] = strand_map[gff_parts[6]]
                if is_gff2:
                    gff_info = _nest_gff2_features(gff_info)
                # features that have parents need to link so we can pick up
                # the relationship
                if gff_info['quals'].has_key('Parent'):
                    final_key = 'child'
                elif gff_info['id']:
                    final_key = 'parent'
                # Handle flat features
                else:
                    final_key = 'feature'
            # otherwise, associate these annotations with the full record
            else:
                final_key = 'annotation'
            return [(final_key, (simplejson.dumps(gff_info) if params.jsonify
                else gff_info))]
    return []

def _gff_line_reduce(map_results, out, params):
    """Reduce part of Map-Reduce; combines results of parsed features.
    """
    final_items = dict()
    for gff_type, final_val in map_results:
        send_val = (simplejson.loads(final_val) if (params.jsonify and 
                    gff_type not in ['directive']) else final_val)
        try:
            final_items[gff_type].append(send_val)
        except KeyError:
            final_items[gff_type] = [send_val]
    for key, vals in final_items.items():
        out.add(key, (simplejson.dumps(vals) if params.jsonify else vals))

class _MultiIDRemapper:
    """Provide an ID remapping for cases where a parent has a non-unique ID.

    Real life GFF3 cases have non-unique ID attributes, which we fix here
    by using the unique sequence region to assign children to the right
    parent.
    """
    def __init__(self, base_id, all_parents):
        self._base_id = base_id
        self._parents = all_parents

    def remap_id(self, feature_dict):
        rstart, rend = feature_dict['location']
        for index, parent in enumerate(self._parents):
            pstart, pend = parent['location']
            if rstart >= pstart and rend <= pend:
                return ("%s_%s" % (self._base_id, index + 1) if index > 0
                        else self._base_id)
        raise ValueError("Did not find remapped ID location: %s, %s, %s" % (
                self._base_id, [p['location'] for p in self._parents],
                feature_dict['location']))

class GFFMapReduceFeatureAdder:
    """Move through a GFF file, adding new features to SeqRecord objects.
    """
    def __init__(self, base_dict=None, line_adjust_fn=None,
            disco_host=None, create_missing=True):
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
        if base_dict is None:
            base_dict = dict()
        self.base = base_dict
        self._create_missing = create_missing
        self._disco_host = disco_host
        self._map_fn = _gff_line_map
        self._reduce_fn = _gff_line_reduce
        self._line_adjust_fn = line_adjust_fn
        # details on what we can filter items with
        self._filter_info = dict(gff_id = [0], gff_source_type = [1, 2],
                gff_type = [2])
    
    def add_features(self, gff_files, limit_info=None):
        """Add all features from a GFF file to our base dictionary.
        """
        [b for b in self.add_features_gen(gff_files, limit_info)]

    def add_features_gen(self, gff_files, limit_info=None, target_lines=None):
        """Generator add_features version which can be used for partial parses.
        """
        # gracefully handle a single file passed
        if isinstance(gff_files, str):
            gff_files = [gff_files]
        # turn all limit information into tuples for identical comparisons
        final_limit_info = {}
        if limit_info:
            for key, values in limit_info.items():
                final_limit_info[key] = [(v,) if isinstance(v, str) 
                        else tuple(v) for v in values]
        if self._disco_host:
            assert target_lines is None, "Cannot split parallelized jobs"
            assert self._line_adjust_fn is None, \
                    "Cannot adjust lines on parallelized jobs"
            results = self._disco_process(gff_files, final_limit_info)
            self._results_to_features(results)
        else:
            for results in self._std_process(gff_files, final_limit_info,
                    target_lines):
                self._results_to_features(results)
                yield

    def _results_to_features(self, results):
        """Add parsed dictionaries of results to Biopython SeqFeatures.
        """
        self._add_annotations(results.get('annotation', []))
        [self._add_toplevel_feature(f) for f in results.get('feature', [])]
        self._add_parent_child_features(results.get('parent', []),
                results.get('child', []))
        self._add_seqs(results.get('fasta', []))
        self._add_directives(results.get('directive', []))

    def _add_directives(self, directives):
        """Handle any directives or meta-data in the GFF file.

        Relevant items are added as annotation meta-data to each record.
        """
        dir_keyvals = collections.defaultdict(list)
        for directive in directives:
            parts = directive.split()
            if len(parts) > 1:
                key = parts[0]
                val = (parts[1] if len(parts) == 2 else tuple(parts[1:]))
                dir_keyvals[key].append(val)
        for key, vals in dir_keyvals.items():
            for rec in self.base.values():
                self._add_ann_to_rec(rec, key, vals)

    def _add_seqs(self, recs):
        """Add sequence information contained in the GFF3 to records.
        """
        for rec in recs:
            if self.base.has_key(rec.id):
                self.base[rec.id].seq = rec.seq
            else:
                self.base[rec.id] = rec
    
    def _add_parent_child_features(self, parents, children):
        """Add nested features with parent child relationships.
        """
        multi_remap = self._identify_dup_ids(parents)
        # add children features
        children_prep = collections.defaultdict(list)
        for child_dict in children:
            child_feature = self._get_feature(child_dict)
            for pindex, pid in enumerate(child_feature.qualifiers['Parent']):
                if multi_remap.has_key(pid):
                    pid = multi_remap[pid].remap_id(child_dict)
                    child_feature.qualifiers['Parent'][pindex] = pid
                children_prep[pid].append((child_dict['rec_id'],
                    child_feature))
        children = dict(children_prep)
        # add children to parents that exist
        for cur_parent_dict in parents:
            cur_id = cur_parent_dict['id']
            if multi_remap.has_key(cur_id):
                cur_parent_dict['id'] = multi_remap[cur_id].remap_id(
                        cur_parent_dict)
            cur_parent = self._add_toplevel_feature(cur_parent_dict)
            cur_parent, children = self._add_children_to_parent(cur_parent,
                    children)
        # create parents for children without them (GFF2 or split/bad files)
        while len(children) > 0:
            parent_id, cur_children = children.items()[0]
            # one child, do not nest it
            if len(cur_children) == 1:
                rec_id, child = cur_children[0]
                loc = (child.location.nofuzzy_start, child.location.nofuzzy_end)
                rec = self._get_rec(dict(rec_id=rec_id, location=loc))
                rec.features.append(child)
                del children[parent_id]
            else:
                cur_parent = self._add_missing_parent(parent_id, cur_children)
                cur_parent, children = self._add_children_to_parent(cur_parent,
                        children)

    def _identify_dup_ids(self, parents):
        """Identify duplicated ID attributes in potential nested parents.

        According to the GFF3 spec ID attributes are supposed to be unique
        for a file, but this is not always true in practice. This looks
        for duplicates, and provides unique IDs sorted by locations.
        """
        multi_ids = collections.defaultdict(list)
        for parent in parents:
            multi_ids[parent['id']].append(parent)
        multi_ids = [(mid, parents) for (mid, parents) in multi_ids.items()
                if len(parents) > 1]
        multi_remap = dict()
        for mid, parents in multi_ids:
            multi_remap[mid] = _MultiIDRemapper(mid, parents)
        return multi_remap

    def _add_children_to_parent(self, cur_parent, children):
        """Recursively add children to parent features.
        """
        if children.has_key(cur_parent.id):
            cur_children = children[cur_parent.id]
            for rec_id, cur_child in cur_children:
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
            rec = self._get_rec(ann)
            for key, vals in ann['quals'].items():
                self._add_ann_to_rec(rec, key, vals)

    def _add_ann_to_rec(self, rec, key, vals):
        """Add a key/value annotation to the given SeqRecord.
        """
        if rec.annotations.has_key(key):
            try:
                rec.annotations[key].extend(vals)
            except AttributeError:
                rec.annotations[key] = [rec.annotations[key]] + vals
        else:
            rec.annotations[key] = vals

    def _get_rec(self, base_dict):
        """Retrieve a record to add features to.
        """
        max_loc = base_dict.get('location', (0, 1))[1]
        try:
            cur_rec = self.base[base_dict['rec_id']]
            # update generated unknown sequences with the expected maximum length
            if isinstance(cur_rec.seq, UnknownSeq):
                cur_rec.seq._length = max([max_loc, cur_rec.seq._length])
            return cur_rec
        except KeyError:
            if self._create_missing:
                new_rec = SeqRecord(UnknownSeq(max_loc), base_dict['rec_id'])
                self.base[base_dict['rec_id']] = new_rec
                return new_rec
            else:
                raise

    def _add_missing_parent(self, parent_id, cur_children):
        """Add a new feature that is missing from the GFF file.
        """
        base_rec_id = list(set(c[0] for c in cur_children))
        assert len(base_rec_id) == 1
        feature_dict = dict(id=parent_id, strand=None,
                type="inferred_parent", quals=dict(), rec_id=base_rec_id[0])
        coords = [(c.location.nofuzzy_start, c.location.nofuzzy_end) 
                for r, c in cur_children]
        feature_dict["location"] = (min([c[0] for c in coords]),
                max([c[1] for c in coords]))
        return self._add_toplevel_feature(feature_dict)

    def _add_toplevel_feature(self, feature_dict):
        """Add a toplevel non-nested feature to the appropriate record.
        """
        new_feature = self._get_feature(feature_dict)
        rec = self._get_rec(feature_dict)
        rec.features.append(new_feature)
        return new_feature

    def _get_feature(self, feature_dict):
        """Retrieve a Biopython feature from our dictionary representation.
        """
        location = FeatureLocation(*feature_dict['location'])
        new_feature = SeqFeature(location, feature_dict['type'],
                id=feature_dict['id'], strand=feature_dict['strand'])
        new_feature.qualifiers = feature_dict['quals']
        return new_feature

    def _get_local_params(self, limit_info=None):
        class _LocalParams:
            def __init__(self):
                self.jsonify = False
        params = _LocalParams()
        params.limit_info = limit_info
        params.filter_info = self._filter_info
        return params

    def _std_process(self, gff_files, limit_info, target_lines):
        """Process GFF addition without any parallelization.

        In addition to limit filtering, this accepts a target_lines attribute
        which provides a number of lines to parse before returning results.
        This allows partial parsing of a file to prevent memory issues.
        """
        params = self._get_local_params(limit_info)
        class _LocalOut:
            def __init__(self, smart_breaks=False):
                self._items = dict()
                self._smart_breaks = smart_breaks
                self._missing_keys = collections.defaultdict(int)
                self.can_break = True
                self.num_lines = 0
            def add(self, key, vals):
                if self._smart_breaks:
                    # if we are not GFF2 we expect parents and break
                    # based on not having missing ones
                    if key == 'directive':
                        if vals[0] == '#':
                            self.can_break = True
                    elif not vals[0].get("is_gff2", False):
                        self._update_missing_parents(key, vals)
                        self.can_break = (len(self._missing_keys) == 0)
                    # otherwise, break more simply when we are done with
                    # stretches of child features
                    else:
                        self.can_break = (key != 'child')
                self.num_lines += 1
                try:
                    self._items[key].extend(vals)
                except KeyError:
                    self._items[key] = vals
            def _update_missing_parents(self, key, vals):
                # smart way of deciding if we can break this.
                # if this is too much, can go back to not breaking in the
                # middle of children
                if key in ["child"]:
                    for val in vals:
                        for p_id in val["quals"]["Parent"]:
                            self._missing_keys[p_id] += 1
                for val in vals:
                    try:
                        del self._missing_keys[val["quals"]["ID"][0]]
                    except KeyError:
                        pass
            def has_items(self):
                return len(self._items) > 0
            def get_results(self):
                return self._items
        out_info = _LocalOut(target_lines is not None)
        for gff_file in gff_files:
            in_handle = open(gff_file)
            found_seqs = False
            while 1:
                line = in_handle.readline()
                if not line:
                    break
                results = self._map_fn(line, params)
                if self._line_adjust_fn and results:
                    if results[0][0] not in ['directive']:
                        results = [(results[0][0],
                            self._line_adjust_fn(results[0][1]))]
                self._reduce_fn(results, out_info, params)
                if (target_lines and out_info.num_lines >= target_lines and
                        out_info.can_break):
                    yield out_info.get_results()
                    out_info = _LocalOut(True)
                if (results and results[0][0] == 'directive' and 
                        results[0][1] == 'FASTA'):
                    found_seqs = True
                    break
            if found_seqs:
                fasta_recs = self._parse_fasta(in_handle)
                out_info.add('fasta', fasta_recs)
            in_handle.close()
        if out_info.has_items():
            yield out_info.get_results()

    def _parse_fasta(self, in_handle):
        """Parse FASTA sequence information contained in the GFF3 file.
        """
        return list(SeqIO.parse(in_handle, "fasta"))

    def _disco_process(self, gff_files, limit_info):
        """Process GFF addition, using Disco to parallelize the process.
        """
        # make these imports local; only need them when using disco
        import simplejson
        import disco
        # absolute path names unless they are special disco files 
        full_files = [(os.path.abspath(f) if f.split(":")[0] != "disco" else f)
                for f in gff_files]
        results = disco.job(self._disco_host, name="gff_reader",
                input=full_files,
                params=disco.Params(limit_info=limit_info, jsonify=True,
                    filter_info=self._filter_info),
                required_modules=["simplejson", "collections", "re"],
                map=self._map_fn, reduce=self._reduce_fn)
        processed = dict()
        for out_key, out_val in disco.result_iterator(results):
            processed[out_key] = simplejson.loads(out_val)
        return processed

class GFFExaminer:
    """Provide high level details about a GFF file to refine parsing.

    GFF is a spec and is provided by many different centers. Real life files
    will present the same information in slightly different ways. Becoming
    familiar with the file you are dealing with is the best way to extract the
    information you need. This class provides high level summary details to
    help in learning.
    """
    def __init__(self):
        feature_adder = GFFMapReduceFeatureAdder()
        self._filter_info = feature_adder._filter_info
        self._local_params = feature_adder._get_local_params()
    
    def available_limits(self, gff_file):
        """Return dictionary information on possible limits for this file.

        This returns a nested dictionary with the following structure:
        
        keys -- names of items to filter by
        values -- dictionary with:
            keys -- filter choice
            value -- counts of that filter in this file

        Not a parallelized map-reduce implementation.
        """
        gff_handle = open(gff_file)
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
        gff_handle.close()
        return final_dict

    def parent_child_map(self, gff_file):
        """Provide a mapping of parent to child relationships in the file.

        Returns a dictionary of parent child relationships:

        keys -- tuple of (source, type) for each parent
        values -- tuple of (source, type) as children of that parent
        
        Not a parallelized map-reduce implementation.
        """
        gff_handle = open(gff_file)
        # collect all of the parent and child types mapped to IDs
        parent_sts = dict()
        child_sts = collections.defaultdict(list)
        for line in gff_handle:
            line_type, line_info = _gff_line_map(line, self._local_params)[0]
            if (line_type == 'parent' or (line_type == 'child' and
                line_info['id'])):
                parent_sts[line_info['id']] = (line_info['quals']['source'][0],
                        line_info['type'])
            if line_type == 'child':
                for parent_id in line_info['quals']['Parent']:
                    child_sts[parent_id].append((line_info['quals']['source'][0],
                        line_info['type']))
        #print parent_sts, child_sts
        # generate a dictionary of the unique final type relationships
        pc_map = collections.defaultdict(list)
        for parent_id, parent_type in parent_sts.items():
            for child_type in child_sts[parent_id]:
                pc_map[parent_type].append(child_type)
        pc_final_map = dict()
        for ptype, ctypes in pc_map.items():
            unique_ctypes = list(set(ctypes))
            unique_ctypes.sort()
            pc_final_map[ptype] = unique_ctypes
        gff_handle.close()
        return pc_final_map

class GFFAddingIterator:
    """Iterate over regions of a GFF file, returning features from each region.
    """
    def __init__(self, seed_dict=None, line_adjust_fn=None,
            feature_adder=None):
        """Initialize with a dictionary of SeqRecords to serve as a seed.

        seed_dict -- dictionary which will be used as the base for every
        group of GFF features added
        feature_adder -- An initialized feature adder like
        GFFMapReduceFeatureAdder, which can be individually customized
        """
        if seed_dict is None:
            seed_dict = dict()
        self._seed = seed_dict
        if feature_adder is None:
            feature_adder = GFFMapReduceFeatureAdder(
                    line_adjust_fn=line_adjust_fn)
        self._adder = feature_adder

    def get_features(self, gff_files, limit_info=None, target_lines=None):
        """Retrieve features from a GFF file added to a seed dict in groups.

        gff_files -- list of one or more GFF files from which features will
        be parsed.
        target_lines -- the number of lines to parse. You are not guaranteed to
        get exactly this many lines as it will be broken to try and keep nested
        features together. If None, the entire set of files will be parsed.
        """
        self._adder.base = copy.deepcopy(self._seed)
        for file_break in self._adder.add_features_gen(gff_files,
                limit_info=limit_info, target_lines=target_lines):
            yield self._adder.base
            self._adder.base = copy.deepcopy(self._seed)

    def get_all_features(self, gff_files, limit_info=None):
        """Retrieve all features from a GFF file, ignoring iterators.
        """
        return self.get_features(gff_files, limit_info).next()
