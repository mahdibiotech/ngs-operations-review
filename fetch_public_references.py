#!/usr/bin/env python3
"""Refresh the compact RefSeq FASTA panel from NCBI E-utilities by accession."""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

from pipeline_v5 import read_fasta, sha256

ROOT = Path(__file__).resolve().parent
META = ROOT / "config/public_reference_metadata.json"


def main():
    metadata = json.loads(META.read_text())
    chunks = []
    for item in metadata["references"]:
        dest = ROOT / item["fasta"]
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
            "db": "nuccore", "id": item["accession"], "rettype": "fasta", "retmode": "text"})
        request = urllib.request.Request(url, headers={"User-Agent": "ngs-lot-review-demo/0.5 contact: research-demo"})
        with urllib.request.urlopen(request, timeout=60) as response:
            text = response.read().decode("utf-8")
        records = list(read_fasta_text(text))
        if len(records) != 1 or records[0][0] != item["accession"]:
            raise RuntimeError(f"Unexpected NCBI FASTA record for {item['accession']}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text if text.endswith("\n") else text+"\n")
        item["length"] = len(records[0][1])
        item["sha256"] = sha256(dest)
        chunks.append(dest.read_text().strip())
        print(f"Fetched {item['accession']}: {item['length']} nt SHA256={item['sha256']}")
        time.sleep(0.35)
    combined = ROOT / metadata["panel_fasta"]
    combined.write_text("\n".join(chunks)+"\n")
    metadata["panel_sha256"] = sha256(combined)
    metadata["retrieved_utc"] = time.strftime("%Y-%m-%d", time.gmtime())
    META.write_text(json.dumps(metadata, indent=2)+"\n")
    print(f"Updated {combined.relative_to(ROOT)} SHA256={metadata['panel_sha256']}")


def read_fasta_text(text):
    name = None; seq = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None: yield name, "".join(seq).upper()
            name = line[1:].split()[0]; seq = []
        elif name is not None:
            seq.append(line.strip())
    if name is not None: yield name, "".join(seq).upper()


if __name__ == "__main__":
    main()
