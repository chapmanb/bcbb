"""Output Biopython SeqRecords and SeqFeatures to GFF3 format.

The target format is GFF3, the current GFF standard:
    http://www.sequenceontology.org/gff3.shtml
"""
import urllib

from Bio import SeqIO

class _IdHandler:
    """Generate IDs for GFF3 Parent/Child relationships where they don't exist.
    """
    def __init__(self):
        self._prefix = "biopygen"
        self._counter = 1
        self._seen_ids = []

    def _generate_id(self, quals):
        """Generate a unique ID not present in our existing IDs.
        """
        gen_id = self._get_standard_id(quals)
        if gen_id is None:
            while 1:
                gen_id = "%s%s" % (self._prefix, self._counter)
                if gen_id not in self._seen_ids:
                    break
                self._counter += 1
        return gen_id

    def _get_standard_id(self, quals):
        """Retrieve standardized IDs from other sources like NCBI GenBank.

        This tries to find IDs from known key/values when stored differently
        than GFF3 specifications.
        """
        possible_keys = ["transcript_id", "protein_id"]
        for test_key in possible_keys:
            if quals.has_key(test_key):
                cur_id = quals[test_key]
                if isinstance(cur_id, tuple) or isinstance(cur_id, list):
                    return cur_id[0]
                else:
                    return cur_id
        return None

    def update_quals(self, quals, has_children):
        """Update a set of qualifiers, adding an ID if necessary.
        """
        cur_id = quals.get("ID", None)
        # if we have an ID, record it
        if cur_id:
            if not isinstance(cur_id, list) and not isinstance(cur_id, tuple):
                cur_id = [cur_id]
            for add_id in cur_id:
                self._seen_ids.append(add_id)
        # if we need one and don't have it, create a new one
        elif has_children:
            new_id = self._generate_id(quals)
            self._seen_ids.append(new_id)
            quals["ID"] = [new_id]
        return quals

class GFF3Writer:
    """Write GFF3 files starting with standard Biopython objects.
    """
    def __init__(self):
        pass

    def write(self, recs, out_handle, include_fasta=False):
        """Write the provided records to the given handle in GFF3 format.
        """
        id_handler = _IdHandler()
        self._write_header(out_handle)
        fasta_recs = []
        try:
            recs = iter(recs)
        except TypeError:
            recs = [recs]
        for rec in recs:
            self._write_rec(rec, out_handle)
            self._write_annotations(rec.annotations, rec.id, len(rec.seq), out_handle)
            for sf in rec.features:
                sf = self._clean_feature(sf)
                id_handler = self._write_feature(sf, rec.id, out_handle,
                        id_handler)
            if include_fasta and len(rec.seq) > 0:
                fasta_recs.append(rec)
        if len(fasta_recs) > 0:
            self._write_fasta(fasta_recs, out_handle)

    def _clean_feature(self, feature):
        quals = {}
        for key, val in feature.qualifiers.items():
            if not isinstance(val, (list, tuple)):
                val = [val]
            val = [str(x) for x in val]
            quals[key] = val
        feature.qualifiers = quals
        clean_sub = [self._clean_feature(f) for f in feature.sub_features]
        feature.sub_features = clean_sub
        return feature

    def _write_rec(self, rec, out_handle):
        # if we have a SeqRecord, write out optional directive
        if len(rec.seq) > 0:
            out_handle.write("##sequence-region %s 1 %s\n" % (rec.id, len(rec.seq)))

    def _get_phase(self, feature):
        if feature.qualifiers.has_key("phase"):
            phase = feature.qualifiers["phase"][0]
        elif feature.type == "CDS":
            phase = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
        else:
            phase = "."
        return str(phase)

    def _write_feature(self, feature, rec_id, out_handle, id_handler,
            parent_id=None):
        """Write a feature with location information.
        """
        if feature.strand == 1:
            strand = '+'
        elif feature.strand == -1:
            strand = '-'
        else:
            strand = '.'
        # remove any standard features from the qualifiers
        quals = feature.qualifiers.copy()
        for std_qual in ["source", "score", "phase"]:
            if quals.has_key(std_qual) and len(quals[std_qual]) == 1:
                del quals[std_qual]
        # add a link to a parent identifier if it exists
        if parent_id:
            if not quals.has_key("Parent"):
                quals["Parent"] = []
            quals["Parent"].append(parent_id)
        quals = id_handler.update_quals(quals, len(feature.sub_features) > 0)
        if feature.type:
            ftype = feature.type
        else:
            ftype = "sequence_feature"
        parts = [str(rec_id),
                 feature.qualifiers.get("source", ["feature"])[0],
                 ftype,
                 str(feature.location.nofuzzy_start + 1), # 1-based indexing
                 str(feature.location.nofuzzy_end),
                 feature.qualifiers.get("score", ["."])[0],
                 strand,
                 self._get_phase(feature),
                 self._format_keyvals(quals)]
        out_handle.write("\t".join(parts) + "\n")
        for sub_feature in feature.sub_features:
            id_handler = self._write_feature(sub_feature, rec_id, out_handle,
                    id_handler, quals["ID"][0])
        return id_handler

    def _format_keyvals(self, keyvals):
        format_kvs = []
        for key in sorted(keyvals.iterkeys()):
            values = keyvals[key]
            key = key.strip()
            format_vals = []
            if not isinstance(values, list) or isinstance(values, tuple):
                values = [values]
            for val in values:
                val = urllib.quote(str(val).strip(), safe=":/ ")
                if ((key and val) and val not in format_vals):
                    format_vals.append(val)
            format_kvs.append("%s=%s" % (key, ",".join(format_vals)))
        return ";".join(format_kvs)

    def _write_annotations(self, anns, rec_id, size, out_handle):
        """Add annotations which refer to an entire sequence.
        """
        format_anns = self._format_keyvals(anns)
        if format_anns:
            parts = [rec_id, "annotation", "remark", "1", str(size if size > 1 else 1),
                     ".", ".", ".", format_anns]
            out_handle.write("\t".join(parts) + "\n")

    def _write_header(self, out_handle):
        """Write out standard header directives.
        """
        out_handle.write("##gff-version 3\n")

    def _write_fasta(self, recs, out_handle):
        """Write sequence records using the ##FASTA directive.
        """
        out_handle.write("##FASTA\n")
        SeqIO.write(recs, out_handle, "fasta")

def write(recs, out_handle, include_fasta=False):
    """High level interface to write GFF3 files from SeqRecords and SeqFeatures.

    If include_fasta is True, the GFF3 file will include sequence information
    using the ##FASTA directive.
    """
    writer = GFF3Writer()
    return writer.write(recs, out_handle, include_fasta)
