#!/usr/bin/env python3
"""Generate a deterministic truth-labeled challenge panel from bundled public FASTA references."""
import argparse
import csv
import random
from pathlib import Path

from pipeline_v5 import read_fasta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", type=Path, default=Path("data/v5_challenge"))
    p.add_argument("--reads-per-sample", type=int, default=5000)
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--length", type=int, default=150)
    p.add_argument("--seed", type=int, default=20260927)
    p.add_argument("--fractions", default="0.001,0.005,0.02,0.2")
    args = p.parse_args()
    if args.reads_per_sample < 100 or args.replicates < 1 or args.length < 50:
        p.error("Use >=100 reads/sample, >=1 replicate, and read length >=50")
    sequences = dict(read_fasta(Path("data/public_references/panel.fasta")))
    viral = ["NC_045512.2", "NC_001803.1"]
    host = sequences["NC_012920.1"]
    for acc in viral:
        if args.length > len(sequences[acc]):
            p.error(f"Read length exceeds reference length for {acc}")
    rng = random.Random(args.seed)
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    manifest = out / "manifest.csv"
    rows = []
    fractions = [float(x) for x in args.fractions.split(",")]
    if any(x < 0 or x > 1 for x in fractions):
        p.error("Fractions must lie between 0 and 1")

    def source_read(sequence):
        start = rng.randrange(0, len(sequence) - args.length + 1)
        seq = sequence[start:start + args.length]
        # Low-rate substitutions produce realistic mismatches without implying a validated model.
        chars = list(seq)
        for i, base in enumerate(chars):
            if base in "ACGT" and rng.random() < 0.002:
                chars[i] = rng.choice([b for b in "ACGT" if b != base])
        return "".join(chars)

    for acc in viral:
        for frac in fractions:
            for rep in range(1, args.replicates + 1):
                sid = f"POS_{acc.replace('.', '_')}_F{frac:g}_R{rep}"
                fastq = out / f"{sid}.fastq"
                virus = sequences[acc]
                n_viral = round(args.reads_per_sample * frac)
                with fastq.open("w") as handle:
                    for i in range(args.reads_per_sample):
                        seq = source_read(virus if i < n_viral else host)
                        handle.write(f"@{sid}_{i+1}\n{seq}\n+\n{'I'*len(seq)}\n")
                rows.append({"sample_id": sid, "fastq": str(fastq), "truth": "positive", "expected_taxon": acc})
    for rep in range(1, args.replicates + 1):
        sid = f"NEG_HOST_R{rep}"
        fastq = out / f"{sid}.fastq"
        with fastq.open("w") as handle:
            for i in range(args.reads_per_sample):
                seq = source_read(host)
                handle.write(f"@{sid}_{i+1}\n{seq}\n+\n{'I'*len(seq)}\n")
        rows.append({"sample_id": sid, "fastq": str(fastq), "truth": "negative", "expected_taxon": ""})
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "fastq", "truth", "expected_taxon"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} labeled samples to {manifest}")
    print("Truth is known because reads were simulated from the indicated bundled references; these are not public clinical samples.")


if __name__ == "__main__":
    main()
