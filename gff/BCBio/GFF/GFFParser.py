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
import itertools

# Make defaultdict compatible with versions of python older than 2.4
try:
    collections.defaultdict
except AttributeError:
    import _utils
    collections.defaultdict = _utils.defaultdict

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
            pieces = []
            for p in parts:
                # fix misplaced semi-colons in keys in some GFF2 files
                if p and p[0] == ';':
                    p = p[1:]
                pieces.append(p.strip().split(" "))
            key_vals = [(p[0], " ".join(p[1:])) for p in pieces]
        for item in key_vals:
            # standard in-spec items are key=value
            if len(item) == 2:
                key, val = item
            # out-of-spec files can have just key values. We set an empty value
            # which will be changed to true later to standardize.
            else:
                assert len(item) == 1, item
                key = item[0]
                val = ''
            # remove quotes in GFF2 files
            if (len(val) > 0 and val[0] == '"' and val[-1] == '"'):
                val = val[1:-1] 
            if val:
                quals[key].extend([v for v in val.split(',') if v])
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
        # case for WormBase GFF -- everything labelled as Transcript or CDS
        for flat_name in ["Transcript", "CDS"]:
            if gff_parts["quals"].has_key(flat_name):
                # parent types
                if gff_parts["type"] in [flat_name]:
                    if not gff_parts["id"]:
                        gff_parts["id"] = gff_parts["quals"][flat_name][0]
                        gff_parts["quals"]["ID"] = [gff_parts["id"]]
                # children types
                elif gff_parts["type"] in ["intron", "exon", "three_prime_UTR",
                        "coding_exon", "five_prime_UTR", "CDS", "stop_codon",
                        "start_codon"]:
                    gff_parts["quals"]["Parent"] = gff_parts["quals"][flat_name]
                break

        return gff_parts

    strand_map = {'+' : 1, '-' : -1, '?' : None, None: None}
    line = line.strip()
    if line[:2] == "##":
        return [('directive', line[2:])]
    elif line and line[0] != "#":
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
            assert len(parts) >= 8, line
            # not python2.4 compatible but easier to understand
            #gff_parts = [(None if p == '.' else p) for p in parts]
            gff_parts = []
            for p in parts:
                if p == ".":
                    gff_parts.append(None)
                else:
                    gff_parts.append(p)
            gff_info = dict()
            # collect all of the base qualifiers for this item
            if len(parts) > 8:
                quals, is_gff2 = _split_keyvals(gff_parts[8])
            else:
                quals, is_gff2 = dict(), False
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
                gff_info['strand'] = strand_map.get(gff_parts[6], None)
                if is_gff2:
                    gff_info = _nest_gff2_features(gff_info)
                # features that have parents need to link so we can pick up
                # the relationship
                if gff_info['quals'].has_key('Parent'):
                    # check for self referential parent/child relationships
                    # remove the ID, which is not useful
                    for p in gff_info['quals']['Parent']:
                        if p == gff_info['id']:
                            gff_info['id'] = ''
                            del gff_info['quals']['ID']
                            break
                    final_key = 'child'
                elif gff_info['id']:
                    final_key = 'parent'
                # Handle flat features
                else:
                    final_key = 'feature'
            # otherwise, associate these annotations with the full record
            else:
                final_key = 'annotation'
            if params.jsonify:
                return [(final_key, simplejson.dumps(gff_info))]
            else:
                return [(final_key, gff_info)]
    return []

def _gff_line_reduce(map_results, out, params):
    """Reduce part of Map-Reduce; combines results of parsed features.
    """
    final_items = dict()
    for gff_type, final_val in map_results:
        if params.jsonify and gff_type not in ['directive']:
            final_val = simplejson.loads(final_val)
        try:
            final_items[gff_type].append(final_val)
        except KeyError:
            final_items[gff_type] = [final_val]
    for key, vals in final_items.items():
        if params.jsonify:
            vals = simplejson.dumps(vals)
        out.add(key, vals)

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
                if index > 0:
                    return ("%s_%s" % (self._base_id, index + 1))
                else:
                    return self._base_id
        raise ValueError("Did not find remapped ID location: %s, %s, %s" % (
                self._base_id, [p['location'] for p in self._parents],
                feature_dict['location']))

class _AbstractMapReduceGFF:
    """Base class providing general GFF parsing for local and remote classes.

    This class should be subclassed to provide a concrete class to parse
    GFF under specific conditions. These classes need to implement
    the _gff_process function, which returns a dictionary of SeqRecord
    information.
    """
    def __init__(self, create_missing=True):
        """Initialize GFF parser 

        create_missing - If True, create blank records for GFF ids not in
        the base_dict. If False, an error will be raised.
        """
        self._create_missing = create_missing
        self._map_fn = _gff_line_map
        self._reduce_fn = _gff_line_reduce
        self._examiner = GFFExaminer()

    def _gff_process(self, gff_files, limit_info, target_lines=None):
        raise NotImplementedError("Derived class must define")

    def parse(self, gff_files, base_dict=None, limit_info=None):
        """Parse a GFF file, returning an iterator of SeqRecords.

        limit_info - A dictionary specifying the regions of the GFF file
        which should be extracted. This allows only relevant portions of a file
        to be parsed.
        
        base_dict - A base dictionary of SeqRecord objects which may be
        pre-populated with sequences and other features. The new features from
        the GFF file will be added to this dictionary.
        """
        for rec in self.parse_in_parts(gff_files, base_dict, limit_info):
            yield rec

    def parse_in_parts(self, gff_files, base_dict=None, limit_info=None,
            target_lines=None):
        """Parse a region of a GFF file specified, returning info as generated.

        target_lines -- The number of lines in the file which should be used
        for each partial parse. This should be determined based on available
        memory.
        """
        for results in self.parse_simple(gff_files, limit_info, target_lines):
            if base_dict is None:
                cur_dict = dict()
            else:
                cur_dict = copy.deepcopy(base_dict)
            cur_dict = self._results_to_features(cur_dict, results)
            all_ids = cur_dict.keys()
            all_ids.sort()
            for cur_id in all_ids:
                yield cur_dict[cur_id]

    def parse_simple(self, gff_files, limit_info=None, target_lines=1):
        """Simple parse which does not build or nest features.

        This returns a simple dictionary representation of each line in the
        GFF file.
        """
        # gracefully handle a single file passed
        if not isinstance(gff_files, (list, tuple)):
            gff_files = [gff_files]
        limit_info = self._normalize_limit_info(limit_info)
        for results in self._gff_process(gff_files, limit_info, target_lines):
            yield results
       
    def _normalize_limit_info(self, limit_info):
        """Turn all limit information into tuples for identical comparisons.
        """
        final_limit_info = {}
        if limit_info:
            for key, values in limit_info.items():
                final_limit_info[key] = []
                for v in values:
                    if isinstance(v, str):
                        final_limit_info[key].append((v,))
                    else:
                        final_limit_info[key].append(tuple(v))
        return final_limit_info
    
    def _results_to_features(self, base, results):
        """Add parsed dictionaries of results to Biopython SeqFeatures.
        """
        base = self._add_annotations(base, results.get('annotation', []))
        for feature in results.get('feature', []):
            (_, base) = self._add_toplevel_feature(base, feature)
        base = self._add_parent_child_features(base, results.get('parent', []),
                results.get('child', []))
        base = self._add_seqs(base, results.get('fasta', []))
        base = self._add_directives(base, results.get('directive', []))
        return base

    def _add_directives(self, base, directives):
        """Handle any directives or meta-data in the GFF file.

        Relevant items are added as annotation meta-data to each record.
        """
        dir_keyvals = collections.defaultdict(list)
        for directive in directives:
            parts = directive.split()
            if len(parts) > 1:
                key = parts[0]
                if len(parts) == 2:
                    val = parts[1]
                else:
                    val = tuple(parts[1:])
                dir_keyvals[key].append(val)
        for key, vals in dir_keyvals.items():
            for rec in base.values():
                self._add_ann_to_rec(rec, key, vals)
        return base

    def _add_seqs(self, base, recs):
        """Add sequence information contained in the GFF3 to records.
        """
        for rec in recs:
            if base.has_key(rec.id):
                base[rec.id].seq = rec.seq
            else:
                base[rec.id] = rec
        return base
    
    def _add_parent_child_features(self, base, parents, children):
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
            cur_parent, base = self._add_toplevel_feature(base, cur_parent_dict)
            cur_parent, children = self._add_children_to_parent(cur_parent,
                    children)
        # create parents for children without them (GFF2 or split/bad files)
        while len(children) > 0:
            parent_id, cur_children = itertools.islice(children.items(),
                    1).next()
            # one child, do not nest it
            if len(cur_children) == 1:
                rec_id, child = cur_children[0]
                loc = (child.location.nofuzzy_start, child.location.nofuzzy_end)
                rec, base = self._get_rec(base,
                        dict(rec_id=rec_id, location=loc))
                rec.features.append(child)
                del children[parent_id]
            else:
                cur_parent, base = self._add_missing_parent(base, parent_id,
                        cur_children)
                cur_parent, children = self._add_children_to_parent(cur_parent,
                        children)
        return base

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
                cur_parent.location_operator = "join"
                cur_parent.sub_features.append(cur_child)
            del children[cur_parent.id]
        return cur_parent, children

    def _add_annotations(self, base, anns):
        """Add annotation data from the GFF file to records.
        """
        # add these as a list of annotations, checking not to overwrite
        # current values
        for ann in anns:
            rec, base = self._get_rec(base, ann)
            for key, vals in ann['quals'].items():
                self._add_ann_to_rec(rec, key, vals)
        return base

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

    def _get_rec(self, base, info_dict):
        """Retrieve a record to add features to.
        """
        max_loc = info_dict.get('location', (0, 1))[1]
        try:
            cur_rec = base[info_dict['rec_id']]
            # update generated unknown sequences with the expected maximum length
            if isinstance(cur_rec.seq, UnknownSeq):
                cur_rec.seq._length = max([max_loc, cur_rec.seq._length])
            return cur_rec, base
        except KeyError:
            if self._create_missing:
                new_rec = SeqRecord(UnknownSeq(max_loc), info_dict['rec_id'])
                base[info_dict['rec_id']] = new_rec
                return new_rec, base
            else:
                raise

    def _add_missing_parent(self, base, parent_id, cur_children):
        """Add a new feature that is missing from the GFF file.
        """
        base_rec_id = list(set(c[0] for c in cur_children))
        assert len(base_rec_id) == 1
        feature_dict = dict(id=parent_id, strand=None,
                type="inferred_parent", quals=dict(ID=[parent_id]),
                rec_id=base_rec_id[0])
        coords = [(c.location.nofuzzy_start, c.location.nofuzzy_end) 
                for r, c in cur_children]
        feature_dict["location"] = (min([c[0] for c in coords]),
                max([c[1] for c in coords]))
        return self._add_toplevel_feature(base, feature_dict)

    def _add_toplevel_feature(self, base, feature_dict):
        """Add a toplevel non-nested feature to the appropriate record.
        """
        new_feature = self._get_feature(feature_dict)
        rec, base = self._get_rec(base, feature_dict)
        rec.features.append(new_feature)
        return new_feature, base

    def _get_feature(self, feature_dict):
        """Retrieve a Biopython feature from our dictionary representation.
        """
        location = FeatureLocation(*feature_dict['location'])
        new_feature = SeqFeature(location, feature_dict['type'],
                id=feature_dict['id'], strand=feature_dict['strand'])
        new_feature.qualifiers = feature_dict['quals']
        return new_feature

    def _parse_fasta(self, in_handle):
        """Parse FASTA sequence information contained in the GFF3 file.
        """
        return list(SeqIO.parse(in_handle, "fasta"))

class _GFFParserLocalOut:
    """Provide a collector for local GFF MapReduce file parsing.
    """
    def __init__(self, smart_breaks=False):
        self._items = dict()
        self._smart_breaks = smart_breaks
        self._missing_keys = collections.defaultdict(int)
        self._last_parent = None
        self.can_break = True
        self.num_lines = 0

    def add(self, key, vals):
        if self._smart_breaks:
            # if we are not GFF2 we expect parents and break
            # based on not having missing ones
            if key == 'directive':
                if vals[0] == '#':
                    self.can_break = True
                self._last_parent = None
            elif not vals[0].get("is_gff2", False):
                self._update_missing_parents(key, vals)
                self.can_break = (len(self._missing_keys) == 0)
            # break when we are done with stretches of child features
            elif key != 'child':
                self.can_break = True
                self._last_parent = None
            # break when we have lots of child features in a row
            # and change between parents
            else:
                cur_parent = vals[0]["quals"]["Parent"][0]
                if (self._last_parent):
                    self.can_break = (cur_parent != self._last_parent)
                self._last_parent = cur_parent
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
        self._last_parent = None
        return self._items

class GFFParser(_AbstractMapReduceGFF):
    """Local GFF parser providing standardized parsing of GFF3 and GFF2 files.
    """
    def __init__(self, line_adjust_fn=None, create_missing=True):
        _AbstractMapReduceGFF.__init__(self, create_missing=create_missing)
        self._line_adjust_fn = line_adjust_fn
    
    def _gff_process(self, gff_files, limit_info, target_lines):
        """Process GFF addition without any parallelization.

        In addition to limit filtering, this accepts a target_lines attribute
        which provides a number of lines to parse before returning results.
        This allows partial parsing of a file to prevent memory issues.
        """
        line_gen = self._file_line_generator(gff_files)
        for out in self._lines_to_out_info(line_gen, limit_info, target_lines):
            yield out

    def _file_line_generator(self, gff_files):
        """Generate single lines from a set of GFF files.
        """
        for gff_file in gff_files:
            if hasattr(gff_file, "read"):
                need_close = False
                in_handle = gff_file
            else:
                need_close = True
                in_handle = open(gff_file)
            found_seqs = False
            while 1:
                line = in_handle.readline()
                if not line:
                    break
                yield line
            if need_close:
                in_handle.close()

    def _lines_to_out_info(self, line_iter, limit_info=None,
            target_lines=None):
        """Generate SeqRecord and SeqFeatures from GFF file lines.
        """
        params = self._examiner._get_local_params(limit_info)
        out_info = _GFFParserLocalOut((target_lines is not None and
                target_lines > 1))
        found_seqs = False
        for line in line_iter:
            results = self._map_fn(line, params)
            if self._line_adjust_fn and results:
                if results[0][0] not in ['directive']:
                    results = [(results[0][0],
                        self._line_adjust_fn(results[0][1]))]
            self._reduce_fn(results, out_info, params)
            if (target_lines and out_info.num_lines >= target_lines and
                    out_info.can_break):
                yield out_info.get_results()
                out_info = _GFFParserLocalOut((target_lines is not None and
                        target_lines > 1))
            if (results and results[0][0] == 'directive' and 
                    results[0][1] == 'FASTA'):
                found_seqs = True
                break

        class FakeHandle:
            def __init__(self, line_iter):
                self._iter = line_iter
            def read(self):
                return "".join(l for l in self._iter)
            def readline(self):
                try:
                    return self._iter.next()
                except StopIteration:
                    return ""

        if found_seqs:
            fasta_recs = self._parse_fasta(FakeHandle(line_iter))
            out_info.add('fasta', fasta_recs)
        if out_info.has_items():
            yield out_info.get_results()

class DiscoGFFParser(_AbstractMapReduceGFF):
    """GFF Parser with parallelization through Disco (http://discoproject.org.
    """
    def __init__(self, disco_host, create_missing=True):
        """Initialize parser.
        
        disco_host - Web reference to a Disco host which will be used for
        parallelizing the GFF reading job.
        """
        _AbstractMapReduceGFF.__init__(self, create_missing=create_missing)
        self._disco_host = disco_host

    def _gff_process(self, gff_files, limit_info, target_lines=None):
        """Process GFF addition, using Disco to parallelize the process.
        """
        assert target_lines is None, "Cannot split parallelized jobs"
        # make these imports local; only need them when using disco
        import simplejson
        import disco
        # absolute path names unless they are special disco files 
        full_files = []
        for f in gff_files:
            if f.split(":")[0] != "disco":
                full_files.append(os.path.abspath(f))
            else:
                full_files.append(f)
        results = disco.job(self._disco_host, name="gff_reader",
                input=full_files,
                params=disco.Params(limit_info=limit_info, jsonify=True,
                    filter_info=self._examiner._filter_info),
                required_modules=["simplejson", "collections", "re"],
                map=self._map_fn, reduce=self._reduce_fn)
        processed = dict()
        for out_key, out_val in disco.result_iterator(results):
            processed[out_key] = simplejson.loads(out_val)
        yield processed

def parse(gff_files, base_dict=None, limit_info=None, target_lines=None):
    """High level interface to parse GFF files into SeqRecords and SeqFeatures.
    """
    parser = GFFParser()
    for rec in parser.parse_in_parts(gff_files, base_dict, limit_info,
            target_lines):
        yield rec

def parse_simple(gff_files, limit_info=None):
    """Parse GFF files as line by line dictionary of parts.
    """
    parser = GFFParser()
    for rec in parser.parse_simple(gff_files, limit_info=limit_info):
        yield rec["child"][0]

def _file_or_handle(fn):
    """Decorator to handle either an input handle or a file.
    """
    def _file_or_handle_inside(*args, **kwargs):
        in_file = args[1]
        if hasattr(in_file, "read"):
            need_close = False
            in_handle = in_file
        else:
            need_close = True
            in_handle = open(in_file)
        args = (args[0], in_handle) + args[2:]
        out = fn(*args, **kwargs)
        if need_close:
            in_handle.close()
        return out
    return _file_or_handle_inside

class GFFExaminer:
    """Provide high level details about a GFF file to refine parsing.

    GFF is a spec and is provided by many different centers. Real life files
    will present the same information in slightly different ways. Becoming
    familiar with the file you are dealing with is the best way to extract the
    information you need. This class provides high level summary details to
    help in learning.
    """
    def __init__(self):
        self._filter_info = dict(gff_id = [0], gff_source_type = [1, 2],
                gff_source = [1], gff_type = [2])
    
    def _get_local_params(self, limit_info=None):
        class _LocalParams:
            def __init__(self):
                self.jsonify = False
        params = _LocalParams()
        params.limit_info = limit_info
        params.filter_info = self._filter_info
        return params
    
    @_file_or_handle
    def available_limits(self, gff_handle):
        """Return dictionary information on possible limits for this file.

        This returns a nested dictionary with the following structure:
        
        keys -- names of items to filter by
        values -- dictionary with:
            keys -- filter choice
            value -- counts of that filter in this file

        Not a parallelized map-reduce implementation.
        """
        cur_limits = dict()
        for filter_key in self._filter_info.keys():
            cur_limits[filter_key] = collections.defaultdict(int)
        for line in gff_handle:
            # when we hit FASTA sequences, we are done with annotations
            if line.startswith("##FASTA"):
                break
            # ignore empty and comment lines
            if line.strip() and line.strip()[0] != "#":
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

    @_file_or_handle
    def parent_child_map(self, gff_handle):
        """Provide a mapping of parent to child relationships in the file.

        Returns a dictionary of parent child relationships:

        keys -- tuple of (source, type) for each parent
        values -- tuple of (source, type) as children of that parent
        
        Not a parallelized map-reduce implementation.
        """
        # collect all of the parent and child types mapped to IDs
        parent_sts = dict()
        child_sts = collections.defaultdict(list)
        for line in gff_handle:
            # when we hit FASTA sequences, we are done with annotations
            if line.startswith("##FASTA"):
                break
            if line.strip():
                line_type, line_info = _gff_line_map(line,
                        self._get_local_params())[0]
                if (line_type == 'parent' or (line_type == 'child' and
                        line_info['id'])):
                    parent_sts[line_info['id']] = (
                            line_info['quals']['source'][0], line_info['type'])
                if line_type == 'child':
                    for parent_id in line_info['quals']['Parent']:
                        child_sts[parent_id].append((
                            line_info['quals']['source'][0], line_info['type']))
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
        return pc_final_map
