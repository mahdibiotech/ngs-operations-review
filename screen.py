#!/usr/bin/env python3
"""Synthetic NGS lot review. Illustrative gates; no diagnostic or GMP use."""
import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

VERSION = '2.0.0-demo'
COLUMNS = ('sample_id', 'taxon', 'reference_accession', 'reads', 'unique_regions',
           'covered_bases', 'reference_bases', 'mean_identity', 'host_similarity_pct')
ROLES = {'test', 'positive', 'negative'}
RULES = ('min_reads', 'min_breadth', 'min_identity_pct', 'min_total_reads',
         'min_unique_regions', 'host_similarity_review_pct')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def integer(value, label, minimum=0):
    if isinstance(value, bool):
        raise ValueError(f'{label}: integer required')
    try:
        number = int(value)
        if str(number) != str(value).strip():
            raise ValueError()
    except (TypeError, ValueError):
        raise ValueError(f'{label}: integer required') from None
    if number < minimum:
        raise ValueError(f'{label}: must be >= {minimum}')
    return number


def decimal(value, label, minimum=0, maximum=100):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f'{label}: number required') from None
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise ValueError(f'{label}: expected {minimum}..{maximum}')
    return number


def read_inputs(hits_path, metadata_path, config_path):
    metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
    config = json.loads(config_path.read_text(encoding='utf-8'))
    if not isinstance(metadata, dict) or not isinstance(config, dict):
        raise ValueError('metadata and config must be JSON objects')
    for key in ('run_id', 'platform', 'assay_context', 'reference_snapshot'):
        if not isinstance(metadata.get(key), str) or not metadata[key].strip():
            raise ValueError(f'metadata: nonempty {key} required')
    samples = metadata.get('samples')
    if not isinstance(samples, dict) or not samples:
        raise ValueError('metadata: nonempty samples object required')
    for sid, sample in samples.items():
        if not isinstance(sid, str) or not sid or not isinstance(sample, dict):
            raise ValueError('metadata: invalid sample')
        if sample.get('role') not in ROLES:
            raise ValueError(f'{sid}: role must be test, positive or negative')
        sample['total_reads'] = integer(sample.get('total_reads'), f'{sid}.total_reads')
        if sample['role'] == 'positive' and (not isinstance(sample.get('expected_taxon'), str)
                                             or not sample['expected_taxon'].strip()):
            raise ValueError(f'{sid}: expected_taxon required for positive control')
    if not {'test', 'positive', 'negative'} <= {s['role'] for s in samples.values()}:
        raise ValueError('metadata: test, positive and negative controls required')
    if set(config) != set(RULES):
        raise ValueError(f'config: expected exactly {", ".join(RULES)}')
    for key in ('min_reads', 'min_total_reads', 'min_unique_regions'):
        config[key] = integer(config[key], f'config.{key}', 1)
    config['min_breadth'] = decimal(config['min_breadth'], 'config.min_breadth', 0, 1)
    for key in ('min_identity_pct', 'host_similarity_review_pct'):
        config[key] = decimal(config[key], f'config.{key}')
    with hits_path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)) or set(reader.fieldnames) != set(COLUMNS):
            raise ValueError('hits: TSV header must contain exactly the documented columns')
        rows = list(reader)
    if not rows:
        raise ValueError('hits: at least one evidence row required')
    seen = set()
    parsed = []
    for line, row in enumerate(rows, 2):
        sid = row['sample_id']
        if sid not in samples or not row['taxon'] or not row['reference_accession'] or None in row:
            raise ValueError(f'hits line {line}: unknown sample, missing value or extra column')
        key = (sid, row['reference_accession'])
        if key in seen:
            raise ValueError(f'hits line {line}: duplicate sample/reference pair')
        seen.add(key)
        item = {'sample_id': sid, 'role': samples[sid]['role'],
                'taxon': row['taxon'], 'reference_accession': row['reference_accession']}
        for field in ('reads', 'unique_regions', 'covered_bases'):
            item[field] = integer(row[field], f'hits line {line}.{field}')
        item['reference_bases'] = integer(row['reference_bases'], f'hits line {line}.reference_bases', 1)
        for field in ('mean_identity', 'host_similarity_pct'):
            item[field] = decimal(row[field], f'hits line {line}.{field}')
        if (item['unique_regions'] > item['reads'] or item['reads'] > samples[sid]['total_reads']
                or item['covered_bases'] > item['reference_bases']
                or (item['reads'] == 0 and (item['covered_bases'] or item['unique_regions']))):
            raise ValueError(f'hits line {line}: inconsistent counts or coverage')
        item['breadth'] = round(item['covered_bases'] / item['reference_bases'], 6)
        item['rpm'] = round(item['reads'] * 1_000_000 / samples[sid]['total_reads'], 3) if samples[sid]['total_reads'] else 0
        parsed.append(item)
    return metadata, config, parsed


def run(hits_path, metadata_path, config_path):
    metadata, config, hits = read_inputs(hits_path, metadata_path, config_path)
    samples = metadata['samples']
    for hit in hits:
        hit['gate_pass'] = (hit['reads'] >= config['min_reads']
                            and hit['breadth'] >= config['min_breadth']
                            and hit['mean_identity'] >= config['min_identity_pct'])
    negative_taxa = {h['taxon'] for h in hits if h['role'] == 'negative' and h['reads'] > 0}
    qc = {
        'expected_positive_detected': all(any(h['sample_id'] == sid and h['taxon'] == sample['expected_taxon']
                                              and h['gate_pass'] and h['unique_regions'] >= config['min_unique_regions']
                                              for h in hits)
                                          for sid, sample in samples.items() if sample['role'] == 'positive'),
        'negative_control_clear': not any(h['role'] == 'negative' and h['gate_pass'] for h in hits),
        'minimum_depth_met': all(s['total_reads'] >= config['min_total_reads'] for s in samples.values()),
    }
    status = 'REVIEW_READY' if all(qc.values()) else 'QC_BLOCKED'
    worklist = []
    for hit in hits:
        flags = []
        if hit['role'] == 'test' and hit['gate_pass']:
            if hit['taxon'] in negative_taxa:
                flags.append('NEGATIVE_BACKGROUND')
            if hit['unique_regions'] < config['min_unique_regions']:
                flags.append('LIMITED_REGION_SUPPORT')
            if hit['host_similarity_pct'] >= config['host_similarity_review_pct']:
                flags.append('HOST_SIMILARITY')
            worklist.append({'sample_id': hit['sample_id'], 'taxon': hit['taxon'],
                             'reference_accession': hit['reference_accession'],
                             'reads': hit['reads'], 'rpm': hit['rpm'],
                             'unique_regions': hit['unique_regions'], 'breadth': hit['breadth'],
                             'identity_pct': hit['mean_identity'],
                             'host_similarity_pct': hit['host_similarity_pct'], 'review_flags': flags,
                             'state': 'ON_HOLD_QC' if status == 'QC_BLOCKED' else 'HUMAN_REVIEW'})
        hit['review_flags'] = flags
    worklist.sort(key=lambda h: (h['sample_id'], h['taxon'], h['reference_accession']))
    return {
        'schema_version': '2.0', 'run_id': metadata['run_id'], 'status': status,
        'assay_context': metadata['assay_context'], 'platform': metadata['platform'],
        'qc': qc, 'samples': samples, 'hits': hits, 'worklist': worklist,
        'demo_rules': config,
        'provenance': {'hits_sha256': digest(hits_path), 'metadata_sha256': digest(metadata_path),
                       'config_sha256': digest(config_path), 'reference_snapshot': metadata['reference_snapshot'],
                       'software': VERSION},
        'limitations': ('Synthetic training data and illustrative thresholds. Flags are triage prompts, '
                        'not evidence of contamination, validated assay results, lot release or GMP compliance.'),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hits', type=Path, default=Path('data/demo_run/hits.tsv'))
    parser.add_argument('--metadata', type=Path, default=Path('data/demo_run/metadata.json'))
    parser.add_argument('--config', type=Path, default=Path('config/demo_rules.json'))
    parser.add_argument('--output', type=Path, default=Path('results/report.json'))
    parser.add_argument('--html', type=Path, default=Path('results/report.html'))
    parser.add_argument('--worklist', type=Path, default=Path('results/worklist.tsv'))
    args = parser.parse_args()
    try:
        report = run(args.hits, args.metadata, args.config)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        from render_html import render
        args.html.parent.mkdir(parents=True, exist_ok=True)
        args.html.write_text(render(report), encoding='utf-8')
        args.worklist.parent.mkdir(parents=True, exist_ok=True)
        with args.worklist.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, delimiter='\t', fieldnames=('sample_id', 'taxon', 'reference_accession', 'reads', 'rpm', 'unique_regions', 'breadth', 'identity_pct', 'host_similarity_pct', 'review_flags', 'state'))
            writer.writeheader()
            for item in report['worklist']:
                writer.writerow({**item, 'review_flags': ','.join(item['review_flags'])})
        print(f'{report["status"]}: {len(report["worklist"])} candidate(s) listed → {args.output}')
        return 0 if report['status'] == 'REVIEW_READY' else 2
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f'Input error: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
