"""Output Biopython SeqRecords and SeqFeatures to GFF3 format.

The target format is GFF3, the current GFF standard:
    http://www.sequenceontology.org/gff3.shtml
"""
import urllib

class GFF3Writer:
    """Write GFF3 files starting with standard Biopython objects.
    """
    def __init__(self):
        pass

    def write(self, recs, out_handle):
        """Write the provided records to the given handle in GFF3 format.
        """
        self._write_header(out_handle)
        for rec in recs:
            self._write_annotations(rec.annotations, rec.id, out_handle)
            for sf in rec.features:
                self._write_feature(sf, rec.id, out_handle)

    def _write_feature(self, feature, rec_id, out_handle):
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
        parts = [str(rec_id),
                 feature.qualifiers.get("source", ["feature"])[0],
                 (feature.type if feature.type else "sequence_feature"),
                 str(feature.location.nofuzzy_start),
                 str(feature.location.nofuzzy_end),
                 feature.qualifiers.get("score", ["."])[0],
                 strand,
                 str(feature.qualifiers.get("phase", ["."])[0]),
                 self._format_keyvals(quals)]
        out_handle.write("\t".join(parts) + "\n")
        for sub_feature in feature.sub_features:
            self._write_feature(sub_feature, rec_id, out_handle)

    def _format_keyvals(self, keyvals):
        format_kvs = []
        for key, values in keyvals.items():
            format_vals = [urllib.quote(v) for v in values]
            format_kvs.append("%s=%s" % (key, ",".join(format_vals)))
        return ";".join(format_kvs)

    def _write_annotations(self, anns, rec_id, out_handle):
        """Add annotations which refer to an entire sequence.
        """
        format_anns = self._format_keyvals(anns)
        if format_anns:
            print repr(format_anns)
            parts = [rec_id, "annotation", "remark", ".", ".", ".", ".", ".",
                     format_anns]
            out_handle.write("\t".join(parts) + "\n")

    def _write_header(self, out_handle):
        """Write out standard header directives.
        """
        out_handle.write("##gff-version 3\n")
