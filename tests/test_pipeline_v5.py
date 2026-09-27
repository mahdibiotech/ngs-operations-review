import tempfile
import unittest
from pathlib import Path

from pipeline_v5 import cigar_positions, parse_blast, parse_sam, performance, read_fasta, read_fastq, wilson


class PipelineV5Tests(unittest.TestCase):
    def test_cigar_reference_coordinates(self):
        self.assertEqual(list(cigar_positions(3, "4M2I3M1D2M")), [2,3,4,5,6,7,8,10,11])

    def test_wilson_small_sample_is_wide(self):
        lo, hi = wilson(2, 2)
        self.assertLess(lo, 0.5)
        self.assertGreater(hi, 0.9)

    def test_metrics_unknown_truth_excluded(self):
        m = performance([
            {"truth":"positive", "detected":True}, {"truth":"positive", "detected":False},
            {"truth":"negative", "detected":False}, {"truth":"negative", "detected":True},
            {"truth":"unknown", "detected":True},
        ])
        self.assertEqual((m["tp"],m["fn"],m["tn"],m["fp"]), (1,1,1,1))
        self.assertEqual(m["unknown_truth_samples_excluded"], 1)
        self.assertEqual(m["sensitivity"], 0.5)

    def test_no_positive_denominator_is_null(self):
        m = performance([{"truth":"negative", "detected":False}])
        self.assertIsNone(m["sensitivity"])
        self.assertIsNone(m["sensitivity_95pct_wilson_ci"])

    def test_fastq_reader_checks_shape(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"x.fastq"
            path.write_text("@r1\nACGT\n+\nIIII\n")
            self.assertEqual(list(read_fastq(path))[0][1], "ACGT")

    def test_fastq_reader_rejects_bad_quality_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"bad.fastq"
            path.write_text("@r1\nACGT\n+\nIII\n")
            with self.assertRaises(ValueError): list(read_fastq(path))

    def test_blast_best_hit_and_cutoffs(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"hits.tsv"
            path.write_text("r1\tVIRUS\t99.0\t100\t100\t1e-30\t180\nr1\tHOST\t98.0\t100\t100\t1e-20\t150\nr2\tVIRUS\t85.0\t100\t100\t1e-10\t100\n")
            result = parse_blast(path, {"VIRUS":{},"HOST":{}}, 90, .8)
            self.assertEqual(result["r1"]["ref"], "VIRUS")
            self.assertNotIn("r2", result)

    def test_blast_refseq_pipe_identifier_is_resolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"hits.tsv"
            path.write_text("r1\tref|VIRUS|\t99.0\t100\t100\t1e-30\t180\n")
            result = parse_blast(path, {"VIRUS":{}}, 90, .8)
            self.assertEqual(result["r1"]["ref"], "VIRUS")

    def test_sam_mapq_and_nonprimary_alignments_are_filtered(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"x.sam"
            path.write_text("@SQ\tSN:VIRUS\tLN:100\n"
                            "r1\t0\tVIRUS\t1\t30\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
                            "r2\t0\tVIRUS\t21\t2\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
                            "r3\t256\tVIRUS\t41\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n")
            counts, covered = parse_sam(path, {"VIRUS":{}}, min_mapq=10)
            self.assertEqual(counts["VIRUS"], 1)
            self.assertEqual(len(covered["VIRUS"]), 10)

    def test_empty_metrics_are_not_zero_claims(self):
        m = performance([])
        self.assertIsNone(m["specificity"])
        self.assertEqual(m["negative_samples"], 0)

    def test_public_reference_panel_is_non_synthetic_and_matches_metadata(self):
        root = Path(__file__).resolve().parents[1]
        refs = dict(read_fasta(root/"data/public_references/panel.fasta"))
        self.assertIn("NC_045512.2", refs)
        self.assertIn("NC_001803.1", refs)
        self.assertIn("NC_012920.1", refs)
        self.assertEqual(len(refs["NC_045512.2"]), 29903)


if __name__ == "__main__":
    unittest.main()
