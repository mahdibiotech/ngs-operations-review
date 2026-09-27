import json
import tempfile
import unittest
from pathlib import Path
from screen import run
from render_html import render

ROOT = Path(__file__).parent
HITS = ROOT / 'data/demo_run/hits.tsv'
META = ROOT / 'data/demo_run/metadata.json'
CONFIG = ROOT / 'config/demo_rules.json'
REFERENCES = ROOT / 'config/synthetic_references.json'


class LotReviewTests(unittest.TestCase):
    def test_normal_run_preserves_ambiguous_signal_for_review(self):
        report = run(HITS, META, CONFIG)
        self.assertEqual(report['status'], 'REVIEW_READY')
        self.assertEqual(len(report['worklist']), 2)
        self.assertEqual(report['worklist'][0]['taxon'], 'SyntheticVirus-A')
        ambiguous = report['worklist'][1]
        self.assertEqual(ambiguous['taxon'], 'SyntheticVirus-C')
        self.assertEqual(ambiguous['review_flags'],
                         ['NEGATIVE_BACKGROUND', 'LIMITED_REGION_SUPPORT', 'HOST_SIMILARITY'])
        self.assertEqual(ambiguous['state'], 'HUMAN_REVIEW')

    def test_negative_control_failure_blocks_all_candidate_review(self):
        report = run(ROOT / 'data/scenarios/negative_failure/hits.tsv', META, CONFIG)
        self.assertEqual(report['status'], 'QC_BLOCKED')
        self.assertFalse(report['qc']['negative_control_clear'])
        self.assertTrue(all(row['state'] == 'ON_HOLD_QC' for row in report['worklist']))

    def test_unexpected_positive_taxon_cannot_satisfy_control(self):
        report = run(ROOT / 'data/scenarios/positive_failure/hits.tsv', META, CONFIG)
        self.assertFalse(report['qc']['expected_positive_detected'])
        self.assertEqual(report['status'], 'QC_BLOCKED')

    def test_insufficient_depth_blocks_lot(self):
        report = run(HITS, ROOT / 'data/scenarios/low_depth/metadata.json', CONFIG)
        self.assertFalse(report['qc']['minimum_depth_met'])
        self.assertEqual(report['status'], 'QC_BLOCKED')

    def test_duplicate_reference_and_nonfinite_metrics_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'hits.tsv'
            lines = HITS.read_text().splitlines()
            path.write_text('\n'.join(lines + [lines[1]]) + '\n')
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                run(path, META, CONFIG)
            path.write_text('\n'.join(lines).replace('98.2', 'NaN') + '\n')
            with self.assertRaisesRegex(ValueError, 'mean_identity'):
                run(path, META, CONFIG)

    def test_provenance_changes_when_evidence_changes(self):
        baseline = run(HITS, META, CONFIG)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'hits.tsv'
            path.write_text(HITS.read_text().replace('24\t5\t1800', '25\t5\t1800'))
            modified = run(path, META, CONFIG)
            self.assertNotEqual(baseline['provenance']['hits_sha256'], modified['provenance']['hits_sha256'])
            self.assertEqual(baseline['provenance']['config_sha256'], modified['provenance']['config_sha256'])

    def test_html_escapes_upstream_taxon(self):
        report = run(HITS, META, CONFIG)
        report['worklist'][0]['taxon'] = '<script>alert(1)</script>'
        html = render(report)
        self.assertIn('&lt;script&gt;', html)
        self.assertNotIn('<script>alert(1)</script>', html)

    def test_reference_catalogue_rejects_mismatches_and_tracks_version(self):
        baseline = run(HITS, META, CONFIG)
        self.assertIn('references_sha256', baseline['provenance'])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'references.json'
            original = REFERENCES.read_text(encoding='utf-8')
            path.write_text(original.replace('"length": 10000', '"length": 9999', 1))
            with self.assertRaisesRegex(ValueError, 'taxon/length mismatch'):
                run(HITS, META, CONFIG, path)
            path.write_text(original.replace('synthetic-refset-2026-09', 'other-snapshot'))
            with self.assertRaisesRegex(ValueError, 'snapshot_id'):
                run(HITS, META, CONFIG, path)
            path.write_text(original.replace('"SYN-D1"', '"SYN-D2"'))
            self.assertNotEqual(baseline['provenance']['references_sha256'],
                                run(HITS, META, CONFIG, path)['provenance']['references_sha256'])

    def test_challenge_matrix_passes(self):
        from validate_demo import replay
        result = replay()
        self.assertTrue(result['all_passed'], result['cases'])
        self.assertEqual(len(result['cases']), 8)

    def test_fastq_demo_reads_to_candidate_report(self):
        from fastq_demo import summarize
        rows, metadata, counts, ambiguous = summarize(
            ROOT / 'data/fastq_demo', ROOT / 'data/fastq_demo/references.fasta',
            ROOT / 'config/fastq_demo_references.json', META)
        self.assertEqual(counts, {'TEST001': 4, 'TEST002': 3, 'NEG001': 1, 'POS001': 4})
        self.assertEqual(ambiguous, dict.fromkeys(counts, 0))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            import csv
            with (path / 'hits.tsv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=('sample_id', 'taxon', 'reference_accession',
                    'reads', 'unique_regions', 'covered_bases', 'reference_bases', 'mean_identity',
                    'host_similarity_pct'), delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)
            (path / 'metadata.json').write_text(json.dumps(metadata))
            report = run(path / 'hits.tsv', path / 'metadata.json',
                         ROOT / 'config/fastq_demo_rules.json',
                         ROOT / 'config/fastq_demo_references.json')
            self.assertEqual(report['status'], 'REVIEW_READY')
            self.assertEqual(len(report['worklist']), 2)

    def test_fastq_demo_rejects_malformed_quality_and_unknown_fasta(self):
        from fastq_demo import read_fastq, summarize
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.fastq'
            path.write_text('@r\nACGT\n+\nIII\n')
            with self.assertRaisesRegex(ValueError, 'malformed'):
                list(read_fastq(path))
            path.write_text('>SYN-X1\nACGT\n')
            with self.assertRaisesRegex(ValueError, 'FASTA catalogue'):
                summarize(ROOT / 'data/fastq_demo', path,
                          ROOT / 'config/fastq_demo_references.json', META)


if __name__ == '__main__':
    unittest.main()
