#!/usr/bin/env python3
"""Educational review of synthetic viral alignment hits; not a diagnostic assay."""
import argparse, csv, hashlib, json, sys
from pathlib import Path

REQUIRED = ('sample_id', 'taxon', 'reads', 'covered_bases', 'reference_bases', 'mean_identity')

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def run(hits_path, metadata_path):
    meta = json.loads(metadata_path.read_text())
    if not isinstance(meta, dict) or not isinstance(meta.get('samples'), dict):
        raise ValueError('metadata: samples object required')
    samples = meta['samples']
    if not samples or any(not isinstance(v, dict) or v.get('role') not in ('test', 'negative', 'positive') or not isinstance(v.get('total_reads'), int) or v['total_reads'] < 0 for v in samples.values()):
        raise ValueError('metadata: each sample requires role and nonnegative integer total_reads')
    roles = {sample['role'] for sample in samples.values()}
    if not {'test', 'negative', 'positive'} <= roles:
        raise ValueError('metadata: test, negative and positive sample roles required')
    with hits_path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames is None or set(REQUIRED) - set(reader.fieldnames):
            raise ValueError('hits: missing required columns')
        hits = list(reader)
    results = []
    seen = set()
    for row in hits:
        if row['sample_id'] not in samples or not row['taxon'].strip():
            raise ValueError('unknown sample or empty taxon')
        key = row['sample_id'], row['taxon']
        if key in seen:
            raise ValueError('duplicate sample/taxon hit')
        seen.add(key)
        reads = int(row['reads']); covered = int(row['covered_bases']); length = int(row['reference_bases']); identity = float(row['mean_identity'])
        if reads < 0 or length <= 0 or not 0 <= covered <= length or not 0 <= identity <= 100:
            raise ValueError('invalid hit metrics')
        coverage = covered / length
        # Demonstration-only review thresholds, not assay acceptance criteria.
        review = reads >= 3 and coverage >= .05 and identity >= 90
        results.append({'sample_id': row['sample_id'], 'role': samples[row['sample_id']]['role'], 'taxon': row['taxon'], 'reads': reads, 'breadth': round(coverage, 4), 'identity_pct': identity, 'review_flag': review})
    positive_ok = all(any(r['sample_id'] == sid and r['review_flag'] for r in results) for sid, sample in samples.items() if sample['role'] == 'positive')
    negative_ok = not any(r['role'] == 'negative' and r['review_flag'] for r in results)
    depth_ok = all(s['total_reads'] >= 1000 for s in samples.values())
    qc = {'positive_control_detected': positive_ok, 'negative_control_clear': negative_ok, 'minimum_read_count_met': depth_ok}
    return {'status': 'REVIEW_READY' if all(qc.values()) else 'QC_BLOCKED', 'qc': qc, 'hits': results, 'provenance': {'hits_sha256': sha256(hits_path), 'metadata_sha256': sha256(metadata_path), 'software': 'pathoquest-demo/0.1', 'purpose': 'synthetic training demonstration'}, 'limitations': 'No clinical or GMP interpretation. Review thresholds are illustrative; a flagged taxon is not a confirmed contaminant.'}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hits', type=Path, required=True); p.add_argument('--metadata', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    try:
        result = run(a.hits, a.metadata)
        a.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
        print(result['status'])
        return 0 if result['status'] == 'REVIEW_READY' else 2
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as exc:
        print(f'Input error: {exc}', file=sys.stderr); return 1
if __name__ == '__main__': sys.exit(main())
