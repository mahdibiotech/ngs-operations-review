#!/usr/bin/env python3
"""Exploratory viral NGS workflow: Bowtie2, BLAST+, optional SPAdes assembly."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "pathoquest-demo/0.5.0"
DEFAULTS = {"min_aligned_reads": 10, "min_blast_reads": 5,
            "min_reference_breadth": 0.05, "min_mapq": 10, "min_read_identity": 90.0,
            "min_read_query_coverage": 0.80, "min_contig_identity": 90.0,
            "min_contig_query_coverage": 0.50}
DNA = set("ACGTN")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path: Path):
    name, seq = None, []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            elif name is None:
                raise ValueError(f"FASTA sequence before header in {path}")
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()


def read_fastq(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
            if not seq or not plus or not qual or not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {path}")
            if len(seq.strip()) != len(qual.strip()) or not set(seq.strip().upper()) <= DNA:
                raise ValueError(f"Invalid sequence/quality in {path}: {header.strip()}")
            yield header[1:].strip().split()[0], seq.strip().upper(), qual.strip()


def fasta_from_fastq(fastq: Path, fasta: Path):
    with fasta.open("w") as out:
        count = 0
        for name, seq, _ in read_fastq(fastq):
            out.write(f">{name}\n{seq}\n")
            count += 1
    if count == 0:
        raise ValueError(f"No reads found in {fastq}")
    return count


def cigar_positions(start1: int, cigar: str):
    """Yield reference coordinates covered by M/= /X operations in a SAM CIGAR."""
    pos = start1 - 1
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(length)
        if op in "M=X":
            yield from range(pos, pos + n)
        if op in "MDN=X":
            pos += n


def parse_sam(path: Path, ref_meta: dict, min_mapq: int = 0):
    counts = defaultdict(int)
    covered = defaultdict(set)
    with path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            flag = int(cols[1])
            if flag & 4 or flag & 256 or flag & 2048 or int(cols[4]) < min_mapq:
                continue
            ref = cols[2]
            if ref not in ref_meta:
                continue
            counts[ref] += 1
            covered[ref].update(cigar_positions(int(cols[3]), cols[5]))
    return counts, covered


def parse_blast(path: Path, ref_meta: dict, min_identity: float, min_qcov: float):
    best = {}
    if not path.exists():
        return best
    with path.open() as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 7:
                continue
            qid, sid, pident, length, qlen, evalue, bitscore = cols
            # BLAST's sseqid can be decorated (for example ref|NC_045512.2|).
            # Resolve either the accession.version field or a token in that form.
            if sid not in ref_meta:
                sid = next((part for part in sid.split("|") if part in ref_meta), sid)
            hit = {"ref": sid, "identity": float(pident), "qcov": int(length) / max(1, int(qlen)),
                   "bitscore": float(bitscore), "evalue": evalue}
            if qid not in best or hit["bitscore"] > best[qid]["bitscore"]:
                best[qid] = hit
    return {q: h for q, h in best.items()
            if h["identity"] >= min_identity and h["qcov"] >= min_qcov
            and h["ref"] in ref_meta}


def wilson(successes: int, total: int, z: float = 1.959963984540054):
    if total == 0:
        return None
    p = successes / total
    den = 1 + z * z / total
    center = (p + z * z / (2 * total)) / den
    half = z * math.sqrt((p * (1-p) + z*z/(4*total)) / total) / den
    return [max(0.0, center-half), min(1.0, center+half)]


def performance(samples):
    labelled = [s for s in samples if s.get("truth") in {"positive", "negative"}]
    tp = sum(s["truth"] == "positive" and s.get("expected_detected", s["detected"]) for s in labelled)
    fn = sum(s["truth"] == "positive" and not s.get("expected_detected", s["detected"]) for s in labelled)
    tn = sum(s["truth"] == "negative" and not s["detected"] for s in labelled)
    fp = sum(s["truth"] == "negative" and s["detected"] for s in labelled)
    return {
        "basis": "sample-level comparison against manifest truth labels; exploratory panel only",
        "tp": tp, "fn": fn, "tn": tn, "fp": fp,
        "sensitivity": tp/(tp+fn) if tp+fn else None,
        "sensitivity_95pct_wilson_ci": wilson(tp, tp+fn),
        "specificity": tn/(tn+fp) if tn+fp else None,
        "specificity_95pct_wilson_ci": wilson(tn, tn+fp),
        "positive_samples": tp+fn, "negative_samples": tn+fp,
        "unknown_truth_samples_excluded": sum(s.get("truth") not in {"positive", "negative"} for s in samples),
    }


def run(cmd, log_path: Path, trace: list, cwd=None):
    exe = shutil.which(cmd[0])
    if exe is None:
        raise RuntimeError(f"Required executable not found: {cmd[0]}. Install the environment in environment.yml.")
    command = [exe, *cmd[1:]]
    cp = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(cp.stdout)
    trace.append({"command": command, "cwd": str(cwd) if cwd else None, "exit_code": cp.returncode, "log": str(log_path)})
    if cp.returncode:
        raise RuntimeError(f"Command failed ({cp.returncode}): {' '.join(command)}; see {log_path}")
    return cp.stdout


def tool_version(program: str, args: list[str]):
    exe = shutil.which(program)
    if exe is None:
        raise RuntimeError(f"Required executable not found: {program}. Install the environment in environment.yml.")
    cp = subprocess.run([exe, *args], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return cp.stdout.strip().splitlines()[0] if cp.returncode == 0 and cp.stdout.strip() else "version unavailable"


def render_report(report: dict, output: Path):
    rows = []
    for s in report["samples"]:
        evidence = []
        for k,v in s["taxa"].items():
            if v["candidate"]:
                evidence.append(f"{html.escape(k)} ({v['aligned_reads']} aligned; {v['blast_reads']} BLAST reads; {v['breadth']:.1%} breadth; {len(v['contig_blast_hits'])} contig hits)")
        assembly = f"Yes ({s['assembly']['contig_count']} contigs; {s['assembly']['blast_hits']} BLAST hits)" if s.get("assembly") else "Not run"
        rows.append(f"<tr><td>{html.escape(s['sample_id'])}</td><td>{html.escape(s.get('truth') or 'unknown')}</td><td>{html.escape('; '.join(evidence) or 'No candidate')}</td><td>{'Detected' if s['detected'] else 'Not detected'}</td><td>{assembly}</td></tr>")
    m = report["performance"]
    def metric(v, ci):
        return "not estimable" if v is None else f"{v:.1%} (95% Wilson CI {ci[0]:.1%}–{ci[1]:.1%})"
    body = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Exploratory viral NGS report</title>
<style>body{{font:16px system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#18212f}}.warn{{background:#fff2cf;padding:1rem;border-left:5px solid #d68a00}}table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #ccd2da;padding:.55rem;text-align:left}}code{{overflow-wrap:anywhere}}</style>
<h1>Exploratory viral NGS analysis report</h1><p><b>Software:</b> {VERSION} · <b>Run:</b> {html.escape(report['run_id'])}</p>
<div class="warn"><b>Research demonstration only.</b> This is not a certificate of analysis, clinical report, lot-release decision, validated diagnostic assay, or GMP/BPF record. Public references are a small panel; decision cutoffs are illustrative. Sensitivity/specificity below describe only this labeled challenge panel.</div>
<h2>Sample results</h2><table><thead><tr><th>Sample</th><th>Panel truth</th><th>Viral candidates</th><th>Screen</th><th>Assembly run</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<h2>Challenge-panel performance</h2><p>Sensitivity (expected taxon detected in positive samples): <b>{metric(m['sensitivity'],m['sensitivity_95pct_wilson_ci'])}</b> (n={m['positive_samples']})</p><p>Specificity (no viral candidate in negative samples): <b>{metric(m['specificity'],m['specificity_95pct_wilson_ci'])}</b> (n={m['negative_samples']})</p><p>TP={m['tp']}, FN={m['fn']}, TN={m['tn']}, FP={m['fp']}. Unknown-truth samples excluded: {m['unknown_truth_samples_excluded']}.</p>
<h2>Methods and provenance</h2><p>Reads are aligned with Bowtie 2, compared against a local BLAST+ nucleotide database, and optionally assembled with SPAdes rnaviral mode. Candidates require the configured alignment count, reference breadth and independent read-BLAST support. A contig BLAST hit adds supporting evidence; it does not by itself establish an organism or contamination.</p>
<p>Reference catalogue SHA-256: <code>{report['reference_catalog_sha256']}</code></p><p>Input and output hashes, software versions, command lines and logs are in the JSON report.</p>
<p><b>Interpretation:</b> a detected candidate is a computational finding for human review. A non-detection applies only to the tested sample, reference panel, library and cutoffs; it is not proof of absence.</p></html>'''
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(body)


def load_manifest(path: Path, root: Path):
    rows = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            fq = (root / row["fastq"]).resolve()
            if not fq.is_file():
                raise ValueError(f"FASTQ does not exist: {fq}")
            truth = row.get("truth", "unknown").strip().lower() or "unknown"
            if truth not in {"positive", "negative", "unknown"}:
                raise ValueError(f"Invalid truth label for {row['sample_id']}: {truth}")
            rows.append({"sample_id": row["sample_id"], "fastq": fq, "truth": truth,
                         "expected_taxon": row.get("expected_taxon", "").strip() or None})
    if not rows:
        raise ValueError("Manifest has no sample rows")
    if len({r["sample_id"] for r in rows}) != len(rows):
        raise ValueError("Sample IDs must be unique")
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True, type=Path, help="CSV: sample_id,fastq,truth,expected_taxon")
    p.add_argument("--references", default=Path("data/public_references/panel.fasta"), type=Path)
    p.add_argument("--reference-metadata", default=Path("config/public_reference_metadata.json"), type=Path)
    p.add_argument("--config", default=Path("config/pipeline_v5.json"), type=Path)
    p.add_argument("--outdir", default=Path("results/v5"), type=Path)
    p.add_argument("--assemble", action="store_true", help="Run SPAdes rnaviral assembly per sample")
    p.add_argument("--assemble-sample", action="append", default=[], help="With --assemble, limit assembly to this sample ID; repeatable")
    a = p.parse_args(argv)
    root = Path.cwd().resolve()
    manifest = a.manifest.resolve(); refs = a.references.resolve(); out = a.outdir.resolve()
    config = json.loads(a.config.read_text()); metadata = json.loads(a.reference_metadata.read_text())
    thresholds = {**DEFAULTS, **config.get("thresholds", {})}
    ref_meta = {r["accession"]: r for r in metadata["references"]}
    fasta_records = list(read_fasta(refs))
    lengths = {name: len(seq) for name, seq in fasta_records}
    if set(lengths) != set(ref_meta):
        raise ValueError("FASTA IDs and metadata accessions do not match exactly")
    manifest_rows = load_manifest(manifest, root)
    for row in manifest_rows:
        expected = row["expected_taxon"]
        if expected and (expected not in ref_meta or ref_meta[expected]["role"] != "virus"):
            raise ValueError(f"expected_taxon must be a viral accession in the selected reference panel: {expected}")
        if row["truth"] == "positive" and not expected:
            raise ValueError(f"Positive truth label requires expected_taxon: {row['sample_id']}")
    manifest_ids = {r["sample_id"] for r in manifest_rows}
    unknown_assembly_ids = set(a.assemble_sample) - manifest_ids
    if unknown_assembly_ids:
        raise ValueError(f"Assembly sample ID(s) not in manifest: {', '.join(sorted(unknown_assembly_ids))}")
    if a.assemble_sample and not a.assemble:
        raise ValueError("--assemble-sample requires --assemble")
    out.mkdir(parents=True, exist_ok=True)
    work = out / "work"; work.mkdir(exist_ok=True)
    trace = []
    bt_index = work / "panel"
    run(["bowtie2-build", str(refs), str(bt_index)], out/"logs/bowtie2-build.log", trace)
    blast_prefix = work / "panel_db"
    run(["makeblastdb", "-in", str(refs), "-dbtype", "nucl", "-parse_seqids", "-out", str(blast_prefix)], out/"logs/makeblastdb.log", trace)
    samples = []
    for row in manifest_rows:
        sid = row["sample_id"]
        sdir = work / sid; sdir.mkdir(exist_ok=True)
        sam = sdir / "aligned.sam"
        run(["bowtie2", "--very-sensitive", "-x", str(bt_index), "-U", str(row["fastq"]), "-S", str(sam)], out/f"logs/{sid}.bowtie2.log", trace)
        aligned, covered = parse_sam(sam, ref_meta, thresholds["min_mapq"])
        query = sdir / "reads.fasta"; read_count = fasta_from_fastq(row["fastq"], query)
        blast_path = sdir / "reads.blast.tsv"
        outfmt = "6 qseqid saccver pident length qlen evalue bitscore"
        run(["blastn", "-task", "megablast", "-query", str(query), "-db", str(blast_prefix), "-outfmt", outfmt, "-max_target_seqs", "5", "-out", str(blast_path)], out/f"logs/{sid}.blast_reads.log", trace)
        read_hits = parse_blast(blast_path, ref_meta, thresholds["min_read_identity"], thresholds["min_read_query_coverage"])
        blast_counts = defaultdict(int)
        for hit in read_hits.values(): blast_counts[hit["ref"]] += 1
        taxa = {}
        viral_detected = False
        for ref, meta in ref_meta.items():
            if meta["role"] != "virus": continue
            breadth = len(covered[ref]) / lengths[ref]
            candidate = (aligned[ref] >= thresholds["min_aligned_reads"] and
                         blast_counts[ref] >= thresholds["min_blast_reads"] and
                         breadth >= thresholds["min_reference_breadth"])
            viral_detected |= candidate
            taxa[meta["label"]] = {"accession": ref, "aligned_reads": aligned[ref],
                "blast_reads": blast_counts[ref], "breadth": breadth,
                "candidate": candidate, "contig_blast_hits": []}
        expected_label = ref_meta[row["expected_taxon"]]["label"] if row["expected_taxon"] in ref_meta else None
        expected_detected = bool(expected_label and expected_label in taxa and taxa[expected_label]["candidate"])
        sample_result = {"sample_id": sid, "truth": row["truth"], "expected_taxon": row["expected_taxon"],
            "read_count": read_count, "detected": viral_detected, "expected_detected": expected_detected, "taxa": taxa,
            "input_sha256": sha256(row["fastq"]), "assembly": None}
        should_assemble = a.assemble and (not a.assemble_sample or sid in a.assemble_sample)
        if should_assemble:
            asm = sdir / "assembly"
            run(["spades.py", "--rnaviral", "--s1", str(row["fastq"]), "-o", str(asm)], out/f"logs/{sid}.spades.log", trace)
            contigs = asm / "contigs.fasta"
            if not contigs.exists(): contigs = asm / "scaffolds.fasta"
            if not contigs.exists():
                raise RuntimeError(f"SPAdes completed but no contigs/scaffolds FASTA was found under {asm}")
            contig_hits = sdir / "contigs.blast.tsv"
            run(["blastn", "-task", "megablast", "-query", str(contigs), "-db", str(blast_prefix), "-outfmt", outfmt,
                 "-max_target_seqs", "5", "-out", str(contig_hits)], out/f"logs/{sid}.blast_contigs.log", trace)
            parsed = parse_blast(contig_hits, ref_meta, thresholds["min_contig_identity"], thresholds["min_contig_query_coverage"])
            for hit in parsed.values():
                meta = ref_meta[hit["ref"]]
                if meta["role"] == "virus": taxa[meta["label"]]["contig_blast_hits"].append(hit)
            sample_result["assembly"] = {"contigs_fasta": str(contigs), "contig_count": sum(1 for _ in read_fasta(contigs)), "blast_hits": len(parsed)}
        samples.append(sample_result)
    versions = {"bowtie2": tool_version("bowtie2", ["--version"]), "blastn": tool_version("blastn", ["-version"]),
                "makeblastdb": tool_version("makeblastdb", ["-version"])}
    if a.assemble: versions["spades"] = tool_version("spades.py", ["--version"])
    report = {"software": VERSION, "tool_versions": versions, "run_id": out.name, "status": "EXPLORATORY_REPORT_READY", "reference_snapshot": metadata["snapshot_id"],
        "reference_catalog_sha256": sha256(refs), "reference_metadata_sha256": sha256(a.reference_metadata),
        "parameter_file_sha256": sha256(a.config), "thresholds": thresholds, "samples": samples,
        "performance": performance(samples), "commands": trace,
        "input_manifest_sha256": sha256(manifest), "limitations": ["Small public reference panel; not comprehensive.",
        "Sensitivity/specificity are estimates on manifest-labeled challenge samples only, with broad uncertainty at small n.",
        "This output is not a certificate of analysis, clinical decision, lot release, or validated GMP/BPF record."]}
    (out/"report.json").write_text(json.dumps(report, indent=2) + "\n")
    render_report(report, out/"report.html")
    print(f"EXPLORATORY_REPORT_READY: {sum(s['detected'] for s in samples)} sample(s) with viral candidate(s)")
    print(f"Sensitivity={report['performance']['sensitivity']} (n={report['performance']['positive_samples']}); Specificity={report['performance']['specificity']} (n={report['performance']['negative_samples']})")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
