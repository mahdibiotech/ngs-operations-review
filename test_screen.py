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


if __name__ == '__main__':
    unittest.main()
