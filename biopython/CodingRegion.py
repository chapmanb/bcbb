"""Represent coding regions, initially for defining changes due to SNPs.
"""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Data import CodonTable

class NonCodingRegion:
    """Represent a standard region of a chromosome or segment, without coding.
    """
    def __init__(self, full_seq, region_name):
        self._seq = full_seq
        self._name = region_name
        self.is_rc = False

    def __str__(self):
        return "Non-coding region: %s" % self._name

    def get_ref_name(self):
        return self._name

    def is_coding(self):
        return False

    def get_feature_details(self):
        return ("", "", "")

    def snp_surround(self, targets, bp_surround):
        """Retrieve the region surrounding a set of SNP targets.
        """
        target_positions = [t['pos'] for t in targets]
        r_start = max(min(target_positions) - bp_surround, 0)
        r_end = min(max(target_positions) + bp_surround, len(self._seq))
        seq_region = self._seq[r_start:r_end]
        targets = [self._add_surround_info(t, r_start) for t in targets]
        return seq_region, targets

    def _add_surround_info(self, snp, region_start):
        """Add local mapping of a SNP to our surrounding alignment region.
        """
        snp["surround_pos"] = snp["pos"] - region_start
        return snp

class CodingRegion:
    """Represent a coding region, providing remapping of coordinates.
    """
    def __init__(self, full_seq, coding_db):
        """Initialize with a sequence and coding database object.
        """
        self._coding_db = coding_db
        table_name = self._coding_db.get("table", "Standard")
        self._surround = 250 #bp
        self._gap = '-'
        if isinstance(full_seq, SeqRecord):
            full_seq = str(full_seq.seq)
        elif isinstance(full_seq, Seq):
            full_seq = str(full_seq)
        self.is_rc = (coding_db['strand'] == -1)
        self._cds_table = CodonTable.unambiguous_dna_by_name[table_name]
        self._local_seq, self._remap, self._upstream, self._downstream = \
                self._build_local_seq(full_seq, coding_db['location'],
                                      coding_db['strand'])
        if coding_db["coding"]:
            self._codons, self._aa_seq = self._build_aa(self._local_seq,
                    self._cds_table, table_name)
        else:
            self._codons = []
            self._aa_seq = ""
        if len(self._codons) == 0:
            self._non_coding_region = NonCodingRegion(full_seq,
                    coding_db["ref_name"])
            self.is_rc = False

    def __str__(self):
        return "%s %s\n\tCoding region: %s" % (self._coding_db["_id"],
                self._coding_db["name"], self._coding_db["location"])
    
    def get_feature_details(self):
        loc_string = ";".join(["%s-%s" % (s, e) for (s, e) in
                               self._coding_db["location"]])
        return (self._coding_db["_id"], self._coding_db["name"], loc_string)

    def get_ref_name(self):
        return self._coding_db["ref_name"]

    def _build_aa(self, seq, cds_table, table_name):
        """Translate our coding sequencing into codons and amino acids.
        """
        if len(self._local_seq) % 3 != 0:
            return ([], "")
        assert seq[:3] in cds_table.start_codons
        aa_seq = str(Seq(seq, unambiguous_dna).translate(table=table_name))
        if (aa_seq.count('*') != 1 or aa_seq[-1] != '*'):
            return ([], "")
        codons = [seq[i*3:(i+1)*3] for i in range(len(seq) // 3)]
        return codons, aa_seq

    def _build_local_seq(self, seq, loc, strand):
        """Generate a local coding sequence with a map of original coordinates.
        """
        remap = {}
        seq_parts = []
        cur_pos = 0
        loc.sort()
        # if we are reverse complemented, add the parts backwards for
        # correct remapping
        if strand == -1:
            loc.reverse()
        for start, end in loc:
            seq_parts.append(seq[start:end])
            cur_region = range(start, end)
            for i, region in enumerate(cur_region):
                if strand == -1:
                    remap[region] = len(cur_region) - 1 - i + cur_pos
                else:
                    remap[region] = i + cur_pos
            cur_pos += len(cur_region)
        if strand == -1:
            seq_parts.reverse()
        lseq = "".join(seq_parts)
        
        upstream = seq[max(loc[0][0] - self._surround, 0):loc[0][0]]
        downstream = seq[loc[-1][1]:min(loc[-1][1] + self._surround, len(seq))]
        if strand == -1:
            lseq = str(Seq(lseq, unambiguous_dna).reverse_complement())
            upstream = str(Seq(downstream,
                unambiguous_dna).reverse_complement())
            downstream = str(Seq(upstream,
                unambiguous_dna).reverse_complement())
        return lseq, remap, upstream, downstream

    def snp_surround(self, targets, bp_surround):
        """Retrieve the codons surrounding a list of target regions.

        This also includes an extra 5' or 3' "codon" with sequence if
        the start or end abuts the end of the sequence.
        """
        # handle special case where we are labelled as coding but really
        # a frameshift or some other non-parseable coding region
        if not self.is_coding():
            return self._non_coding_region.snp_surround(targets, bp_surround)

        targets = [self._add_local_info(t) for t in targets]
        target_positions = [t['codon_pos'] for t in targets]
        num_surround = bp_surround // 3
        up_extra, down_extra = (None, None)
        r_start = min(target_positions) - num_surround
        if r_start < 0:
            up_extra = self._upstream[r_start*3:]
            r_start = 0
        r_end = max(target_positions) + num_surround
        if r_end > len(self._codons):
            down_extra = self._downstream[:(r_end - len(self._codons))*3]
            r_end = len(self._codons)
        codons = self._codons[r_start:r_end]
        target_indexes = [t - r_start for t in target_positions]
        if up_extra:
            up_fake_codons = [up_extra[i*3:(i+1)*3] for 
                    i in range(len(up_extra) // 3)]
            codons = up_fake_codons + codons
            target_indexes = [i + len(up_fake_codons) for i in target_indexes]
        if down_extra:
            codons.insert(-1, down_extra)
        targets = [self._add_surround_info(t, target_indexes[i]) for i, t in
                enumerate(targets)]
        return "".join(codons), targets

    def is_coding(self):
        return len(self._codons) > 0

    def _add_surround_info(self, target, index_pos):
        """Convert a codon index position in a surround region into a position.
        """
        target["surround_pos"] = (index_pos * 3) + target["in_codon_pos"]
        return target

    def _add_local_info(self, snp_info):
        """Add info to a SNP about its position within this coding region.
        """
        local_pos = self._remap[snp_info['pos']]
        ori_base = snp_info['ref_base']
        new_base = snp_info['snp_base']
        if self.is_rc:
            ori_base = str(Seq(ori_base, unambiguous_dna).reverse_complement())
            new_base = str(Seq(new_base, unambiguous_dna).reverse_complement())
        snp_info['codon_pos'] = local_pos // 3
        snp_info['in_codon_pos'] = local_pos % 3
        orig_codon = self._codons[snp_info['codon_pos']]
        mod_codon = list(orig_codon)
        # substitution or deletion
        if ori_base != self._gap:
            assert self._local_seq[local_pos] == ori_base, (
                    self._local_seq[local_pos], ori_base)
            assert orig_codon[snp_info['in_codon_pos']] == ori_base
            mod_codon[snp_info['in_codon_pos']] = new_base
        # insertion
        else:
            codon_pos = snp_info['in_codon_pos']
            mod_codon = mod_codon[:codon_pos] + [new_base] + \
                    mod_codon[codon_pos:]
        mod_codon = "".join(mod_codon)
        snp_info['orig_codon'] = orig_codon
        snp_info['new_codon'] = mod_codon
        return snp_info

    def get_aa(self, codon):
        return ("*" if codon in self._cds_table.stop_codons else
                self._cds_table.forward_table[codon])

