import json, tempfile, unittest
from pathlib import Path
from screen import run

ROOT = Path(__file__).parent
class ScreenTests(unittest.TestCase):
    def test_baseline(self):
        output = run(ROOT/'data/hits.tsv', ROOT/'data/metadata.json')
        self.assertEqual(output['status'], 'REVIEW_READY')
        self.assertEqual([r['taxon'] for r in output['hits'] if r['role']=='test' and r['review_flag']], ['SyntheticVirus-A'])
    def test_contaminated_negative_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'hits.tsv'
            path.write_text((ROOT/'data/hits.tsv').read_text().replace('NEG001\tSyntheticVirus-A\t0\t0\t10000\t0', 'NEG001\tSyntheticVirus-A\t4\t600\t10000\t95'))
            self.assertEqual(run(path, ROOT/'data/metadata.json')['status'], 'QC_BLOCKED')
    def test_missing_positive_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'hits.tsv'
            path.write_text((ROOT/'data/hits.tsv').read_text().replace('POS001\tSyntheticVirus-A\t35\t2500\t10000\t99.2', 'POS001\tSyntheticVirus-A\t1\t100\t10000\t99.2'))
            self.assertEqual(run(path, ROOT/'data/metadata.json')['status'], 'QC_BLOCKED')
    def test_invalid_coverage_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'hits.tsv'
            path.write_text((ROOT/'data/hits.tsv').read_text().replace('800\t10000', '11000\t10000'))
            with self.assertRaises(ValueError): run(path, ROOT/'data/metadata.json')
    def test_duplicate_hit_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'hits.tsv'
            lines = (ROOT/'data/hits.tsv').read_text().splitlines()
            path.write_text('\n'.join(lines + [lines[1]]) + '\n')
            with self.assertRaisesRegex(ValueError, 'duplicate'): run(path, ROOT/'data/metadata.json')
    def test_missing_control_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'metadata.json'
            meta = json.loads((ROOT/'data/metadata.json').read_text())
            del meta['samples']['NEG001']
            path.write_text(json.dumps(meta))
            with self.assertRaisesRegex(ValueError, 'roles required'): run(ROOT/'data/hits.tsv', path)
if __name__ == '__main__': unittest.main()
