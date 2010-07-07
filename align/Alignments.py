"""Retrieve regions of alignments and calculate conservation.
"""
import os
import subprocess
import StringIO
import collections

from bx.align import maf
import numpy

from Phast import PhastConsCommandline

class AlignConservationFinder:
    """Provide organism specific statistics on conservation from a alignment.
    """
    def __init__(self, work_dir):
        self._work_dir = work_dir
        self._gap = '-'
        self._base_params = {
            "--target-coverage" : "0.125",
            "--expected-length" : "15",
            }

    def _get_estimate_aligns(self, rec, orgs, chroms, retriever):
        """Retrieve the subset of alignments for use in estimating parameters.
        """
        # 7.8M -> 71m
        # 703M -> 813m
        target_size = 250000 #250kb
        find_size = 1000 # 1kb
        est_aligns = []
        for target_start in range(0, len(rec), target_size):
            aligns = retriever.get_regions(orgs, chroms, target_start, 
                    target_start + find_size, do_slice=False)
            for a in aligns:
                if a not in est_aligns:
                    est_aligns.append(a)
        return est_aligns

    def estimate_phastcons_trees(self, rec, orgs, chroms, retriever, init_model_file):
        """Use a subset of alignments in a sequence to estimate phastCons trees.

        phastCons will estimate conserved and non-conserved tree parameters from
        our full data set. This needs to be done once for each model file, and
        is repeated only when the file is missing.
        """
        base_name = init_model_file.replace("-init.mod", "")
        (base, _) = os.path.splitext(os.path.split(base_name)[-1])
        tree_root = base_name + "-trees"
        cons_tree_file = "%s.cons.mod" % tree_root
        noncons_tree_file = "%s.noncons.mod" % tree_root
        if (not os.path.exists(cons_tree_file) or
                not os.path.exists(noncons_tree_file)):
            aligns = self._get_estimate_aligns(rec, orgs, chroms, retriever)
            estimate_align_file = self._write_maf_file(
                    os.path.join(self._work_dir, "%s-estimate.maf" % base), aligns)
            self._run_phastcons_estimate(estimate_align_file, tree_root,
                    init_model_file)
        return cons_tree_file, noncons_tree_file

    def conservation_stats(self, base_orgs, chroms, align, cons_model,
            noncons_model):
        """Provide conservation statistics for alignments relative to a base.

        base_orgs is a list of synonyms for the organism to be used as a
        base. It is assumed to be present a single time in all alignments.
        """
        base_org, cmp_orgs = self._retrieve_components(base_orgs, chroms, align)
        conserved_file, score_file = self._run_phastcons_conservation(align,
                cons_model, noncons_model)
        c_scores = self._read_conservation_scores(score_file)
        c_regions = self._read_conserved_regions(conserved_file, align)
        assert len(c_scores) == len(base_org.text.replace(self._gap, "")), \
                (len(c_scores), len(base_org.text))
        score_window = 20
        windows = []
        for i in range(max(1, len(c_scores) - score_window)):
            windows.append(numpy.median(c_scores[i:i+score_window]))
        print numpy.median(c_scores), max(windows), c_regions
        for to_remove in [conserved_file, score_file]:
            if os.path.exists(to_remove):
                os.remove(to_remove)
        return numpy.median(c_scores), max(windows), c_regions

    def _read_conserved_regions(self, in_file, align):
        regions = []
        with open(in_file) as in_handle:
            for line in in_handle:
                parts = line.split()
                start = int(parts[1])
                end = int(parts[2])
                slice_c = align.components[0].slice_by_coord(start, end)
                regions.append(slice_c.text.replace(self._gap, '').upper())
        return regions

    def _read_conservation_scores(self, score_file):
        with open(score_file) as in_handle:
            header = in_handle.readline()
            scores = [float(l) for l in in_handle]
        return scores

    def _run_phastcons_estimate(self, align_file, tree_root, init_model_file):
        """Run phastCons, estimating conserved and non-conserved parameters.
        """
        estimate_params = {
            "--msa-format" : "MAF",
            "alignment" : align_file,
            "--estimate-trees" : tree_root,
            "--no-post-probs" : True,
            "models" : init_model_file
            }
        estimate_params.update(self._base_params)
        estimate_cl = PhastConsCommandline(**estimate_params)
        p = subprocess.Popen(str(estimate_cl).split(), stderr=subprocess.PIPE)
        p.wait()
        error_str = p.stderr.read()
        if error_str.find("ERROR") >= 0:
            print error_str
            raise ValueError("Problem running phastCons")

    def _run_phastcons_conservation(self, align, cons_model, non_cons_model):
        """Run phastCons calculating conservation scores
        """
        base_name = os.path.join(self._work_dir, "phastcons-run")
        conserved_file = base_name + "-conserved.bed"
        score_file = base_name + "-scores.bed"
        align_file = self._write_maf_file(base_name + ".maf", [align])
        calcluate_params = {
            "--msa-format" : "MAF",
            "alignment" : align_file,
            "--most-conserved" : conserved_file,
            "models" : "%s,%s" % (cons_model, non_cons_model),
            }
        calcluate_params.update(self._base_params)
       
        calcluate_cl = PhastConsCommandline(**calcluate_params)
        #print calcluate_cl
        error_out = StringIO.StringIO()
        with open(score_file, "w") as score_handle:
            p = subprocess.Popen(str(calcluate_cl).split(), stdout=score_handle,
                    stderr=subprocess.PIPE)
            p.wait()
        error_str = p.stderr.read()
        if error_str.find("ERROR") >= 0:
            print error_str
            raise ValueError("Problem running phastCons")
        for to_remove in [align_file]:
            if os.path.exists(to_remove):
                os.remove(to_remove)
        return conserved_file, score_file

    def _write_maf_file(self, out_file, aligns):
        with open(out_file, "w") as out_handle:
            writer = maf.Writer(out_handle)
            for align in aligns:
                writer.write(align)
        return out_file

    def _retrieve_components(self, base_orgs, chroms, align):
        base_org = None
        cmp_orgs = []
        for c in align.components:
            c_org, c_chrom = c.src.split(".")[:2]
            if c_chrom in chroms:
                if c_org in base_orgs:
                    assert base_org is None
                    base_org = c
                else:
                    cmp_orgs.append(c)
        assert base_org is not None
        return base_org, cmp_orgs

    def _percent_match(self, subs_table):
        """Calculate simple statistics on matches and gaps.
        """
        gap = 0
        match = 0
        mismatch = 0
        for (base_one, base_two), count in subs_table.items():
            if base_one == base_two:
                if base_one != self._gap:
                    match += 1
            elif base_one == self._gap or base_two == self._gap:
                gap += 1
            else:
                mismatch += 1
        total = float(gap + match + mismatch)
        return (float(match) / total, float(gap) / total)

    def _get_substitutions(self, base_org, cmp_org):
        """Retrieve statistics for substitions between the two organisms.
        """
        base_seq = base_org.text.upper()
        cmp_seq = cmp_org.text.upper()
        assert len(base_seq) == len(cmp_seq)
        subs_table = collections.defaultdict(int)
        for i, b_base in enumerate(base_seq):
            c_base = cmp_seq[i]
            subs_table[(b_base, c_base)] += 1
        return dict(subs_table)

class AlignRetriever:
    """Retrieve multiple alignments corresponding to a chromosome segment.

    This uses a bx-python alignment index to retrieve alignment regions
    corresponding to provided coordinates.
    """
    def __init__(self, index):
        self._index = index
        self._gap = '-'

    def get_regions(self, orgs, chroms, start, end, do_slice=True):
        final_aligns = []
        for org in orgs:
            for chrom in chroms:
                region_id = "%s.%s" % (org, chrom)
                aligns = self._index.get(region_id, start, end)
                for align in aligns:
                    region_start, region_end, region_ori = \
                            self._find_region_start(region_id, align)
                    # if our reference strand is reversed, re-orient this
                    # relative to the forward strand
                    if region_ori == "-":
                        align = align.reverse_complement()
                    gap_remap = self._get_gap_remap(region_id, align)
                    # if we are only partially overlapping, we start with our region
                    cur_start = max(start, region_start)
                    rel_start = cur_start - region_start
                    cur_end = min(end, region_end - 1)
                    rel_end = cur_end - region_start
                    if do_slice:
                        align = align.slice(gap_remap[rel_start],
                                gap_remap[rel_end])
                    falign = self._organize_components(orgs, chroms, align)
                    if falign.text_size > 0:
                        final_aligns.append(falign)
        return final_aligns
    
    def _organize_components(self, base_orgs, chroms, align):
        """Separate our base organism and comparisons from the alignment.
        """
        base_org = None
        cmp_orgs = []
        for c in align.components:
            c_org, c_chrom = c.src.split(".")[:2]
            if c_chrom in chroms:
                if c_org in base_orgs:
                    assert base_org is None
                    base_org = c
                else:
                    # remove any all N alignments, which are useless
                    bases = list(set(c.text.replace(self._gap, '').upper()))
                    if not(len(bases) == 1 and bases[0] == 'N'):
                        cmp_orgs.append(c)
        assert base_org is not None
        align.components = [base_org] + cmp_orgs
        align.remove_all_gap_columns()
        return align

    def _get_gap_remap(self, region_id, align):
        for c in align.components:
            if c.src == region_id:
                return self._remap_gaps(0, c.text)

    def _remap_gaps(self, start_index, orig_seq, gap="-"):
        """Generate a dictionary mapping original coordinates to a gapped alignment.
        """
        pos_remap = dict()
        orig_i = start_index
        for align_i, base in enumerate(orig_seq):
            align_i += start_index
            pos_remap[orig_i] = align_i
            if base != "-":
                orig_i += 1
        return pos_remap

    def _find_region_start(self, region_id, align):
        """Find the chromosomal start of the current region in this alignment.

        Assumes 1 region per set of alignments for the organism.
        """
        for c in align.components:
            if c.src == region_id:
                return (c.forward_strand_start, c.forward_strand_end,
                        c.strand)
        raise ValueError("Did not find region in alignment: %s" % region_id)
