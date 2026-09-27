#!/usr/bin/env python3
"""Small synthetic FASTQ to evidence demonstration; exact matches only."""
import argparse
import csv
import json
from pathlib import Path

from screen import digest, run

ROOT = Path(__file__).resolve().parent
ALPHABET = set('ACGT')
FIELDS = ('sample_id', 'taxon', 'reference_accession', 'reads', 'unique_regions',
          'covered_bases', 'reference_bases', 'mean_identity', 'host_similarity_pct')


def read_fasta(path):
    records = {}
    name, parts = None, []
    for line in path.read_text(encoding='utf-8').splitlines() + ['>END']:
        if line.startswith('>'):
            if name is not None:
                seq = ''.join(parts).upper()
                if not seq or set(seq) - ALPHABET:
                    raise ValueError(f'FASTA: invalid sequence {name}')
                records[name] = seq
            name, parts = line[1:].split()[0], []
            if name == 'END':
                break
            if name in records:
                raise ValueError(f'FASTA: repeated accession {name}')
        elif not name or not line.strip():
            raise ValueError('FASTA: invalid record')
        else:
            parts.append(line.strip())
    if not records:
        raise ValueError('FASTA: no references')
    return records


def read_fastq(path):
    with path.open(encoding='utf-8') as handle:
        while header := handle.readline():
            seq, plus, qual = (handle.readline().rstrip('\n\r') for _ in range(3))
            seq = seq.upper().rstrip('\n\r')
            if (not header.startswith('@') or plus != '+' or not seq or set(seq) - ALPHABET
                    or len(qual) != len(seq) or any(not 33 <= ord(c) <= 126 for c in qual)):
                raise ValueError(f'FASTQ: malformed read in {path}')
            yield seq


def revcomp(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def summarize(fastq_dir, fasta_path, catalogue_path, metadata_path, region_size=20):
    reference_sequences = read_fasta(fasta_path)
    catalogue = json.loads(catalogue_path.read_text(encoding='utf-8'))
    metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
    if catalogue['snapshot_id'] != metadata['reference_snapshot']:
        raise ValueError('FASTA catalogue snapshot mismatch')
    for accession, entry in catalogue['references'].items():
        if accession not in reference_sequences or len(reference_sequences[accession]) != entry['length']:
            raise ValueError(f'FASTA catalogue length/accession mismatch: {accession}')
    if set(reference_sequences) != set(catalogue['references']):
        raise ValueError('FASTA: unregistered accession')
    evidence, counts, ambiguous = [], {}, {}
    for sid, sample in metadata['samples'].items():
        path = fastq_dir / f'{sid}.fastq'
        if not path.is_file():
            raise ValueError(f'FASTQ missing for {sid}')
        hits = {accession: {'reads': 0, 'regions': set(), 'positions': set()}
                for accession in reference_sequences}
        total, ambiguous[sid] = 0, 0
        for seq in read_fastq(path):
            total += 1
            matches = []
            for accession, reference in reference_sequences.items():
                positions = {i for query in {seq, revcomp(seq)}
                             for i in range(len(reference) - len(query) + 1)
                             if reference.startswith(query, i)}
                if positions:
                    matches.append((accession, positions))
            if len(matches) != 1 or len(matches[0][1]) != 1:
                ambiguous[sid] += bool(matches)
                continue
            accession, positions = matches[0]
            hit = hits[accession]
            hit['reads'] += 1
            pos = next(iter(positions))
            hit['positions'].update(range(pos, pos + len(seq)))
            hit['regions'].add(pos // region_size)
        counts[sid] = total
        sample['total_reads'] = total
        for accession, hit in hits.items():
            if hit['reads'] or sample['role'] != 'test':
                evidence.append({'sample_id': sid, 'taxon': catalogue['references'][accession]['taxon'],
                                 'reference_accession': accession, 'reads': hit['reads'],
                                 'unique_regions': min(hit['reads'], len(hit['regions'])),
                                 'covered_bases': len(hit['positions']),
                                 'reference_bases': len(reference_sequences[accession]),
                                 'mean_identity': 100 if hit['reads'] else 0,
                                 'host_similarity_pct': 'NA'})
    return evidence, metadata, counts, ambiguous


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastq-dir', type=Path, default=ROOT / 'data/fastq_demo')
    parser.add_argument('--fasta', type=Path, default=ROOT / 'data/fastq_demo/references.fasta')
    parser.add_argument('--catalogue', type=Path, default=ROOT / 'config/fastq_demo_references.json')
    parser.add_argument('--metadata', type=Path, default=ROOT / 'data/demo_run/metadata.json')
    parser.add_argument('--rules', type=Path, default=ROOT / 'config/fastq_demo_rules.json')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results/fastq_demo')
    args = parser.parse_args()
    try:
        rows, metadata, counts, ambiguous = summarize(args.fastq_dir, args.fasta, args.catalogue, args.metadata)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        hits = args.output_dir / 'hits.tsv'
        with hits.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        meta = args.output_dir / 'metadata.json'
        meta.write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
        result = run(hits, meta, args.rules, args.catalogue)
        result['upstream'] = {'method': 'exact substring, unique accession only; reverse complement included',
                              'fasta_sha256': digest(args.fasta),
                              'fastq_sha256': {sid: digest(args.fastq_dir / f'{sid}.fastq') for sid in counts},
                              'ambiguous_reads': ambiguous, 'read_counts': counts,
                              'note': 'mean_identity=100 for exact matches; host_similarity_pct is NA (not measured).'}
        (args.output_dir / 'report.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
        from render_html import render
        (args.output_dir / 'report.html').write_text(render(result), encoding='utf-8')
        print(f'{result["status"]}: {len(result["worklist"])} candidates; {args.output_dir / "report.html"}')
        return 0 if result['status'] == 'REVIEW_READY' else 2
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        parser.exit(1, f'Input error: {exc}\n')


if __name__ == '__main__':
    raise SystemExit(main())
