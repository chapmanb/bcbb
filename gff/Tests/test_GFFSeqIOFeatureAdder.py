"""Test decoration of existing SeqRecords with GFF through a SeqIO interface.
"""
import sys
import os
import unittest
import pprint
import StringIO

from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio.GFF import (GFF3Writer, GFFExaminer, GFFParser, DiscoGFFParser)

class MapReduceGFFTest(unittest.TestCase):
    """Tests GFF parsing using a map-reduce framework for parallelization.
    """
    def setUp(self):
        self._test_dir = os.path.join(os.path.dirname(__file__), "GFF")
        self._test_gff_file = os.path.join(self._test_dir,
                "c_elegans_WS199_shortened_gff.txt")
        self._disco_host = "http://localhost:7000"
    
    def t_local_map_reduce(self):
        """General map reduce framework without parallelization.
        """
        cds_limit_info = dict(
                gff_type = ["gene", "mRNA", "CDS"],
                gff_id = ['I']
                )
        rec_dict = SeqIO.to_dict(GFF.parse(self._test_gff_file,
            limit_info=cds_limit_info))
        test_rec = rec_dict['I']
        assert len(test_rec.features) == 32

    def t_disco_map_reduce(self):
        """Map reduce framework parallelized using disco.
        """
        # this needs to be more generalized but fails okay with no disco
        try:
            import disco
            import simplejson
        except ImportError:
            print "Skipping -- disco and json not found"
            return
        cds_limit_info = dict(
                gff_source_type = [('Non_coding_transcript', 'gene'),
                             ('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        parser = DiscoGFFParser(disco_host=self._disco_host)
        rec_dict = SeqIO.to_dict(parser.parse(self._test_gff_file,
            limit_info=cds_limit_info))
        final_rec = rec_dict['I']
        # second gene feature is multi-parent
        assert len(final_rec.features) == 2 # two gene feature

class GFF3Test(unittest.TestCase):
    """Real live GFF3 tests from WormBase and NCBI.

    Uses GFF3 data from:

    ftp://ftp.wormbase.org/pub/wormbase/genomes/c_elegans/
    genome_feature_tables/GFF3/
    ftp://ftp.wormbase.org/pub/wormbase/genomes/c_elegans/sequences/dna/

    and from NCBI.
    """
    def setUp(self):
        self._test_dir = os.path.join(os.path.dirname(__file__), "GFF")
        self._test_seq_file = os.path.join(self._test_dir,
                "c_elegans_WS199_dna_shortened.fa")
        self._test_gff_file = os.path.join(self._test_dir,
                "c_elegans_WS199_shortened_gff.txt")
        self._test_gff_ann_file = os.path.join(self._test_dir,
                "c_elegans_WS199_ann_gff.txt")
        self._full_dir = "/usr/home/chapmanb/mgh/ruvkun_rnai/wormbase/" + \
                "data_files_WS198"
        self._test_ncbi = os.path.join(self._test_dir,
                "ncbi_gff3.txt")

    def not_t_full_celegans(self):
        """Test the full C elegans chromosome and GFF files.

        This is used to test GFF on large files and is not run as a standard
        test. You will need to download the files and adjust the paths
        to run this.
        """
        # read the sequence information
        seq_file = os.path.join(self._full_dir, "c_elegans.WS199.dna.fa")
        gff_file = os.path.join(self._full_dir, "c_elegans.WS199.gff3")
        seq_handle = open(seq_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        seq_handle.close()
        #with open(gff_file) as gff_handle:
        #    possible_limits = feature_adder.available_limits(gff_handle)
        #    pprint.pprint(possible_limits)
        rnai_types = [('Orfeome', 'PCR_product'),
                    ('GenePair_STS', 'PCR_product'),
                    ('Promoterome', 'PCR_product')]
        gene_types = [('Non_coding_transcript', 'gene'),
                      ('Coding_transcript', 'gene'),
                      ('Coding_transcript', 'mRNA'),
                      ('Coding_transcript', 'CDS')]
        limit_info = dict(gff_source_type = rnai_types + gene_types)
        for rec in GFF.parse(gff_file, seq_dict, limit_info=limit_info):
            pass

    def _get_seq_dict(self):
        """Internal reusable function to get the sequence dictionary.
        """
        seq_handle = open(self._test_seq_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        seq_handle.close()
        return seq_dict
    
    def t_possible_limits(self):
        """Calculate possible queries to limit a GFF file.
        """
        gff_examiner = GFFExaminer()
        possible_limits = gff_examiner.available_limits(self._test_gff_file)
        print
        pprint.pprint(possible_limits)

    def t_parent_child(self):
        """Summarize parent-child relationships in a GFF file.
        """
        gff_examiner = GFFExaminer()
        pc_map = gff_examiner.parent_child_map(self._test_gff_file)
        print
        pprint.pprint(pc_map)

    def t_flat_features(self):
        """Check addition of flat non-nested features to multiple records.
        """
        seq_dict = self._get_seq_dict()
        pcr_limit_info = dict(
            gff_source_type = [('Orfeome', 'PCR_product'),
                         ('GenePair_STS', 'PCR_product'),
                         ('Promoterome', 'PCR_product')]
            )
        parser = GFFParser()
        rec_dict = SeqIO.to_dict(parser.parse(self._test_gff_file, seq_dict,
            limit_info=pcr_limit_info))
        assert len(rec_dict['I'].features) == 4
        assert len(rec_dict['X'].features) == 5

    def t_nested_features(self):
        """Check three-deep nesting of features with gene, mRNA and CDS.
        """
        seq_dict = self._get_seq_dict()
        cds_limit_info = dict(
                gff_source_type = [('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        parser = GFFParser()
        rec_dict = SeqIO.to_dict(parser.parse(self._test_gff_file, seq_dict,
            limit_info=cds_limit_info))
        final_rec = rec_dict['I']
        # first gene feature is plain
        assert len(final_rec.features) == 2 # two gene feature
        assert len(final_rec.features[0].sub_features) == 1 # one transcript
        # 15 final CDS regions
        assert len(final_rec.features[0].sub_features[0].sub_features) == 15

    def t_nested_multiparent_features(self):
        """Verify correct nesting of features with multiple parents.
        """
        seq_dict = self._get_seq_dict()
        cds_limit_info = dict(
                gff_source_type = [('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        parser = GFFParser()
        rec_dict = SeqIO.to_dict(parser.parse(self._test_gff_file, seq_dict,
            limit_info=cds_limit_info))
        final_rec = rec_dict['I']
        # second gene feature is multi-parent
        assert len(final_rec.features) == 2 # two gene feature
        cur_subs = final_rec.features[1].sub_features
        assert len(cur_subs) == 3 # three transcripts
        # the first and second transcript have the same CDSs
        assert len(cur_subs[0].sub_features) == 6
        assert len(cur_subs[1].sub_features) == 6
        assert cur_subs[0].sub_features[0] is cur_subs[1].sub_features[0]

    def t_no_dict_error(self):
        """Ensure an error is raised when no dictionary to map to is present.
        """
        parser = GFFParser(create_missing=False)
        try:
            for rec in parser.parse(self._test_gff_file):
                pass
            # no error -- problem
            raise AssertionError('Did not complain with missing dictionary')
        except KeyError:
            pass

    def t_unknown_seq(self):
        """Prepare unknown base sequences with the correct length.
        """
        rec_dict = SeqIO.to_dict(GFF.parse(self._test_gff_file))
        assert len(rec_dict["I"].seq) == 12766937
        assert len(rec_dict["X"].seq) == 17718531

    def t_gff_annotations(self):
        """Check GFF annotations placed on an entire sequence.
        """
        parser = GFFParser()
        rec_dict = SeqIO.to_dict(parser.parse(self._test_gff_ann_file))
        final_rec = rec_dict['I']
        assert len(final_rec.annotations.keys()) == 2
        assert final_rec.annotations['source'] == ['Expr_profile']
        assert final_rec.annotations['expr_profile'] == ['B0019.1']
    
    def t_gff3_iterator(self):
        """Iterated parsing in GFF3 files with nested features.
        """
        parser = GFFParser()
        recs = [r for r in parser.parse_in_parts(self._test_gff_file,
            target_lines=70)]
        # should be one big set because we don't have a good place to split
        assert len(recs) == 6
        assert len(recs[0].features) == 59
    
    def t_gff3_iterator_limit(self):
        """Iterated interface using a limit query on GFF3 files.
        """
        cds_limit_info = dict(
                gff_source_type = [('Coding_transcript', 'gene'),
                             ('Coding_transcript', 'mRNA'),
                             ('Coding_transcript', 'CDS')],
                gff_id = ['I']
                )
        parser = GFFParser()
        rec_dict = SeqIO.to_dict(parser.parse(self._test_gff_file,
            limit_info=cds_limit_info))
        assert len(rec_dict) == 1
        tfeature = rec_dict["I"].features[0].sub_features[0]
        for sub_test in tfeature.sub_features:
            assert sub_test.type == "CDS", sub_test

    def t_gff3_noval_attrib(self):
        """Parse GFF3 file from NCBI with a key/value pair with no value.
        """
        parser = GFFParser()
        rec_dict = SeqIO.to_dict(parser.parse(self._test_ncbi))
        assert len(rec_dict) == 1
        t_feature = rec_dict.values()[0].features[0]
        assert t_feature.qualifiers["pseudo"] == ["true"]

    def t_gff3_multiple_ids(self):
        """Deal with GFF3 with non-unique ID attributes, using NCBI example.
        """
        parser = GFFParser()
        rec_dict = SeqIO.to_dict(parser.parse(self._test_ncbi))
        assert len(rec_dict) == 1
        t_features = rec_dict.values()[0].features[1:]
        # 4 feature sets, same ID, different positions, different attributes
        assert len(t_features) == 4
        for f in t_features:
            assert len(f.sub_features) == 3

    def t_simple_parsing(self):
        """Parse GFF into a simple line by line dictionary without nesting.
        """
        parser = GFFParser()
        num_lines = 0
        for line_info in parser.parse_simple(self._test_gff_file):
            num_lines += 1
        assert num_lines == 177, num_lines
        line_info = line_info['child'][0]
        assert line_info['quals']['confirmed_est'] == \
                ['yk1055g06.5', 'OSTF085G5_1']
        assert line_info['location'] == [4582718, 4583189]

    def t_extra_comma(self):
        """Correctly handle GFF3 files with extra trailing commas.
        """
        tfile = os.path.join(self._test_dir, "mouse_extra_comma.gff3")
        in_handle = open(tfile)
        for rec in GFF.parse(in_handle):
            pass
        in_handle.close()
        tested = False
        for sub_top in rec.features[0].sub_features:
            for sub in sub_top.sub_features:
                if sub.qualifiers.get("Name", "") == ["CDS:NC_000083.5:LOC100040603"]:
                    tested = True
                    assert len(sub.qualifiers["Parent"]) == 1
        assert tested, "Did not find sub-feature to test"

    def t_novalue_key(self):
        """Handle GFF3 files with keys and no values.
        """
        tfile = os.path.join(self._test_dir, "glimmer_nokeyval.gff3")
        rec = GFF.parse(tfile).next()
        f1, f2 = rec.features
        assert f1.qualifiers['ID'] == ['GL0000006']
        assert len(f1.sub_features) == 2
        assert f1.sub_features[0].qualifiers["Lack 3'-end"] == ["true"]
        assert not f1.sub_features[0].qualifiers.has_key("ID")
        assert f2.qualifiers["Complete"] == ["true"]

class SolidGFFTester(unittest.TestCase):
    """Test reading output from SOLiD analysis, as GFF3.

    See more details on SOLiD GFF here:

    http://solidsoftwaretools.com/gf/project/matogff/
    """
    def setUp(self):
        self._test_dir = os.path.join(os.path.dirname(__file__), "GFF")
        self._test_gff_file = os.path.join(self._test_dir,
                "F3-unique-3.v2.gff")

    def t_basic_solid_parse(self):
        """Basic parsing of SOLiD GFF results files.
        """
        parser = GFFParser()
        rec_dict = SeqIO.to_dict(parser.parse(self._test_gff_file))
        test_feature = rec_dict['3_341_424_F3'].features[0]
        assert test_feature.location.nofuzzy_start == 102716
        assert test_feature.location.nofuzzy_end == 102736
        assert len(test_feature.qualifiers) == 7
        assert test_feature.qualifiers['score'] == ['10.6']
        assert test_feature.qualifiers['source'] == ['solid']
        assert test_feature.strand == -1
        assert test_feature.type == 'read'
        assert test_feature.qualifiers['g'] == ['T2203031313223113212']
        assert len(test_feature.qualifiers['q']) == 20
    
    def t_solid_iterator(self):
        """Iterated parsing in a flat file without nested features.
        """
        parser = GFFParser()
        feature_sizes = []
        for rec in parser.parse_in_parts(self._test_gff_file,
                target_lines=5):
            feature_sizes.append(len(rec.features))
        assert len(feature_sizes) == 112
        assert max(feature_sizes) == 1

    def t_line_adjust(self):
        """Adjust lines during parsing to fix potential GFF problems.
        """
        def adjust_fn(results):
            rec_index = results['quals']['i'][0]
            read_name = results['rec_id']
            results['quals']['read_name'] = [read_name]
            results['rec_id'] = rec_index
            return results
        parser = GFFParser(line_adjust_fn=adjust_fn)
        recs = [r for r in parser.parse(self._test_gff_file)]
        assert len(recs) == 1
        work_rec = recs[0]
        assert work_rec.id == '1'
        assert len(work_rec.features) == 112
        assert work_rec.features[0].qualifiers['read_name'] == \
                ['3_336_815_F3']

class GFF2Tester(unittest.TestCase):
    """Parse GFF2 and GTF files, building features.
    """
    def setUp(self):
        self._test_dir = os.path.join(os.path.dirname(__file__), "GFF")
        self._ensembl_file = os.path.join(self._test_dir, "ensembl_gtf.txt")
        self._wormbase_file = os.path.join(self._test_dir, "wormbase_gff2.txt")
        self._jgi_file = os.path.join(self._test_dir, "jgi_gff2.txt")
        self._wb_alt_file = os.path.join(self._test_dir,
                "wormbase_gff2_alt.txt")

    def t_basic_attributes(self):
        """Parse out basic attributes of GFF2 from Ensembl GTF.
        """
        limit_info = dict(
                gff_source_type = [('snoRNA', 'exon')]
                )
        rec_dict = SeqIO.to_dict(GFF.parse(self._ensembl_file,
            limit_info=limit_info))
        work_rec = rec_dict['I']
        assert len(work_rec.features) == 1
        test_feature = work_rec.features[0]
        qual_keys = test_feature.qualifiers.keys()
        qual_keys.sort()
        assert qual_keys == ['Parent', 'exon_number', 'gene_id', 'gene_name',
                'source', 'transcript_id', 'transcript_name']
        assert test_feature.qualifiers['source'] == ['snoRNA']
        assert test_feature.qualifiers['transcript_name'] == ['NR_001477.2']
        assert test_feature.qualifiers['exon_number'] == ['1']

    def t_tricky_semicolons(self):
        """Parsing of tricky semi-colon positions in WormBase GFF2.
        """
        limit_info = dict(
                gff_source_type = [('Genomic_canonical', 'region')]
                )
        rec_dict = SeqIO.to_dict(GFF.parse(self._wormbase_file,
            limit_info=limit_info))
        work_rec = rec_dict['I']
        assert len(work_rec.features) == 1
        test_feature = work_rec.features[0]
        assert test_feature.qualifiers['Note'] == \
          ['Clone cTel33B; Genbank AC199162', 'Clone cTel33B; Genbank AC199162']

    def t_jgi_gff(self):
        """Parsing of JGI formatted GFF2, nested using transcriptId and proteinID
        """
        rec_dict = SeqIO.to_dict(GFF.parse(self._jgi_file))
        tfeature = rec_dict['chr_1'].features[0]
        assert tfeature.location.nofuzzy_start == 37060
        assert tfeature.location.nofuzzy_end == 38216
        assert tfeature.type == 'inferred_parent'
        assert len(tfeature.sub_features) == 6
        sfeature = tfeature.sub_features[1]
        assert sfeature.qualifiers['proteinId'] == ['873']
        assert sfeature.qualifiers['phase'] == ['0']

    def t_ensembl_nested_features(self):
        """Test nesting of features with GFF2 files using transcript_id.
        """
        rec_dict = SeqIO.to_dict(GFF.parse(self._ensembl_file))
        assert len(rec_dict["I"].features) == 2
        t_feature = rec_dict["I"].features[0]
        assert len(t_feature.sub_features) == 32

    def t_wormbase_nested_features(self):
        """Test nesting of features with GFF2 files using Transcript only.
        """
        rec_dict = SeqIO.to_dict(GFF.parse(self._wormbase_file))
        assert len(rec_dict) == 3
        parent_features = [f for f in rec_dict["I"].features if f.type ==
                "Transcript"]
        assert len(parent_features) == 1
        inferred_features = [f for f in rec_dict["I"].features if f.type ==
                "inferred_parent"]
        assert len(inferred_features) == 0
        tfeature = parent_features[0]
        assert tfeature.qualifiers["WormPep"][0] == "WP:CE40797"
        assert len(tfeature.sub_features) == 46

    def t_wb_cds_nested_features(self):
        """Nesting of GFF2 features with a flat CDS key value pair.
        """
        rec_dict = SeqIO.to_dict(GFF.parse(self._wb_alt_file))
        assert len(rec_dict) == 2
        features = rec_dict.values()[1].features
        assert len(features) == 1
        tfeature = features[0]
        assert tfeature.id == "cr01.sctg102.wum.2.1"
        assert len(tfeature.sub_features) == 7

    def t_gff2_iteration(self):
        """Test iterated features with GFF2 files, breaking without parents.
        """
        recs = []
        for rec in GFF.parse(self._wormbase_file, target_lines=15):
            recs.append(rec)
        assert len(recs) == 4
        assert recs[0].features[0].type == 'region'
        assert recs[0].features[1].type == 'SAGE_tag'
        assert len(recs[0].features[2].sub_features) == 29

class DirectivesTest(unittest.TestCase):
    """Tests for parsing directives and other meta-data.
    """
    def setUp(self):
        self._test_dir = os.path.join(os.path.dirname(__file__), "GFF")
        self._gff_file = os.path.join(self._test_dir, "hybrid1.gff3")

    def t_basic_directives(self):
        """Parse out top level meta-data supplied in a GFF3 file.
        """

        recs = SeqIO.to_dict(GFF.parse(self._gff_file))
        anns = recs['chr17'].annotations
        assert anns['gff-version'] == ['3']
        assert anns['attribute-ontology'] == ['baz']
        assert anns['feature-ontology'] == ['bar']
        assert anns['source-ontology'] == ['boo']
        assert anns['sequence-region'] == [('foo', '1', '100'), ('chr17',
            '62467934', '62469545')]

    def t_fasta_directive(self):
        """Parse FASTA sequence information contained in a GFF3 file.
        """
        recs = SeqIO.to_dict(GFF.parse(self._gff_file))
        assert len(recs) == 1
        test_rec = recs['chr17']
        assert str(test_rec.seq) == "GATTACAGATTACA"
    
    def t_examiner_with_fasta(self):
        """Perform high level examination of files with FASTA directives.
        """
        examiner = GFFExaminer()
        pc_map = examiner.parent_child_map(self._gff_file)
        assert pc_map[('UCSC', 'mRNA')] == [('UCSC', 'CDS')]
        limits = examiner.available_limits(self._gff_file)
        assert limits['gff_id'].keys()[0][0] == 'chr17'
        assert sorted(limits['gff_source_type'].keys()) == \
                [('UCSC', 'CDS'), ('UCSC', 'mRNA')]

class OutputTest(unittest.TestCase):
    """Tests to write SeqFeatures to GFF3 output format.
    """
    def setUp(self):
        self._test_dir = os.path.join(os.path.dirname(__file__), "GFF")
        self._test_seq_file = os.path.join(self._test_dir,
                "c_elegans_WS199_dna_shortened.fa")
        self._test_gff_file = os.path.join(self._test_dir,
                "c_elegans_WS199_shortened_gff.txt")
        self._test_gff_ann_file = os.path.join(self._test_dir,
                "c_elegans_WS199_ann_gff.txt")
        self._wormbase_file = os.path.join(self._test_dir, "wormbase_gff2.txt")

    def t_gff3_to_gff3(self):
        """Read in and write out GFF3 without any loss of information.
        """
        recs = SeqIO.to_dict(GFF.parse(self._test_gff_file))
        out_handle = StringIO.StringIO()
        GFF.write(recs.values(), out_handle)
        wrote_handle = StringIO.StringIO(out_handle.getvalue())
        recs_two = SeqIO.to_dict(GFF.parse(wrote_handle))

        orig_rec = recs.values()[0]
        re_rec = recs.values()[0]
        assert len(orig_rec.features) == len(re_rec.features)
        for i, orig_f in enumerate(orig_rec.features):
            assert str(orig_f) == str(re_rec.features[i])

    def t_gff2_to_gff3(self):
        """Read in GFF2 and write out as GFF3.
        """
        recs = SeqIO.to_dict(GFF.parse(self._wormbase_file))
        out_handle = StringIO.StringIO()
        GFF.write(recs.values(), out_handle)
        wrote_handle = StringIO.StringIO(out_handle.getvalue())
        # check some tricky lines in the GFF2 file
        checks = 0
        for line in wrote_handle:
            if line.find("Interpolated_map_position") >= 0:
                checks += 1
                assert line.find("RFLP=No") > 0
            if line.find("Gene=WBGene00000138") > 0:
                checks += 1
                assert line.find("ID=B0019.1") > 0
            if line.find("translated_nucleotide_match\t12762127") > 0:
                checks += 1
                assert line.find("Note=MSP%3AFADFSPLDVSDVNFATDDLAK") > 0
        assert checks == 3, "Missing check line"

    def t_write_from_recs(self):
        """Write out GFF3 from SeqRecord inputs.
        """
        seq = Seq("GATCGATCGATCGATCGATC")
        rec = SeqRecord(seq, "ID1")
        qualifiers = {"source": "prediction", "score": 10.0, "other": ["Some", "annotations"],
                      "ID": "gene1"}
        sub_qualifiers = {"source": "prediction"}
        top_feature = SeqFeature(FeatureLocation(0, 20), type="gene", strand=1,
                                                          qualifiers=qualifiers)
        top_feature.sub_features = [SeqFeature(FeatureLocation(0, 5), type="exon", strand=1,
                                               qualifiers=sub_qualifiers),
                                    SeqFeature(FeatureLocation(15, 20), type="exon", strand=1,
                                               qualifiers=sub_qualifiers)]
        rec.features = [top_feature]
        out_handle = StringIO.StringIO()
        GFF.write([rec], out_handle)
        wrote_info = out_handle.getvalue().split("\n")
        assert wrote_info[0] == "##gff-version 3"
        assert wrote_info[1] == "##sequence-region ID1 1 20"
        assert wrote_info[2].split("\t") == ['ID1', 'prediction', 'gene', '1',
                                             '20', '10.0', '+', '.',
                                             'other=Some,annotations;ID=gene1']
        assert wrote_info[3].split("\t") == ['ID1', 'prediction', 'exon', '1', '5',
                                             '.', '+', '.', 'Parent=gene1']

    def t_write_fasta(self):
        """Include FASTA records in GFF output.
        """
        seq = Seq("GATCGATCGATCGATCGATC")
        rec = SeqRecord(seq, "ID1")
        qualifiers = {"source": "prediction", "score": 10.0, "other": ["Some", "annotations"],
                      "ID": "gene1"}
        rec.features = [SeqFeature(FeatureLocation(0, 20), type="gene", strand=1,
                                   qualifiers=qualifiers)]
        out_handle = StringIO.StringIO()
        GFF.write([rec], out_handle, include_fasta=True)
        wrote_info = out_handle.getvalue().split("\n")
        fasta_parts = wrote_info[3:]
        assert fasta_parts[0] == "##FASTA"
        assert fasta_parts[1] == ">ID1 <unknown description>"
        assert fasta_parts[2] == str(seq)

    def t_write_seqrecord(self):
        """Write single SeqRecords.
        """
        seq = Seq("GATCGATCGATCGATCGATC")
        rec = SeqRecord(seq, "ID1")
        qualifiers = {"source": "prediction", "score": 10.0, "other": ["Some", "annotations"],
                      "ID": "gene1"}
        rec.features = [SeqFeature(FeatureLocation(0, 20), type="gene", strand=1,
                                   qualifiers=qualifiers)]
        out_handle = StringIO.StringIO()
        GFF.write([rec], out_handle, include_fasta=True)
        wrote_info = out_handle.getvalue().split("\n")
        gff_line = wrote_info[2]
        assert gff_line.split("\t")[0] == "ID1"

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    test_suite = unittest.TestSuite()
    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [GFF3Test, MapReduceGFFTest, SolidGFFTester, GFF2Tester,
             DirectivesTest, OutputTest]
    #tests = [GFF3Test]
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)
    return test_suite

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
