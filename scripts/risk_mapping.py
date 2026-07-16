#!/usr/bin/env python3
"""Map RISK / GSE57945 samples to diagnosis + SRA run accession.

GEO series matrix gives GSM, diagnosis, and the SRA experiment (SRX), aligned by column.
ENA (reliable) resolves the study's runs (run <-> experiment). Output: risk_mapping.tsv
(run, experiment_srx, gsm, diagnosis) + diagnosis counts.
"""
import gzip
import re
import urllib.request
from collections import Counter

MTX = "gse.mtx.gz"  # GSE57945_series_matrix.txt.gz, pre-downloaded


def fetch(url):
    return urllib.request.urlopen(url, timeout=180).read().decode("utf-8", "replace")


# 1. GEO matrix -> GSM, diagnosis, SRX (aligned columns)
txt = gzip.open(MTX, "rt", encoding="utf-8", errors="replace").read()
gsm = diag = srx = None
for line in txt.splitlines():
    if line.startswith("!Sample_geo_accession"):
        gsm = [c.strip().strip('"') for c in line.split("\t")[1:]]
    elif line.startswith("!Sample_characteristics_ch1"):
        v = [c.strip().strip('"') for c in line.split("\t")[1:]]
        if v and v[0].lower().startswith("diagnosis:"):
            diag = [x.split(":", 1)[1].strip() if ":" in x else "" for x in v]
    elif line.startswith("!Sample_relation") and "SRA" in line:
        v = [c.strip().strip('"') for c in line.split("\t")[1:]]
        srx = [(re.search(r"SRX\d+", x).group(0) if re.search(r"SRX\d+", x) else "")
               for x in v]

srx2info = {s: (g, d) for g, s, d in zip(gsm, srx, diag) if s}

# 2. ENA: study accession from one SRX, then all runs for the study
first = next(s for s in srx if s)
lines = fetch("https://www.ebi.ac.uk/ena/portal/api/filereport"
              f"?accession={first}&result=read_run"
              "&fields=study_accession,run_accession,experiment_accession"
              "&format=tsv").splitlines()
h = lines[0].split("\t")
study = lines[1].split("\t")[h.index("study_accession")].strip()

tsv = fetch("https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={study}&result=read_run"
            "&fields=run_accession,experiment_accession&format=tsv&limit=0").splitlines()
h = tsv[0].split("\t")
ri, ei = h.index("run_accession"), h.index("experiment_accession")

out = []
for line in tsv[1:]:
    f = line.split("\t")
    run, exp = f[ri], f[ei]
    g, d = srx2info.get(exp, ("", ""))
    out.append((run, exp, g, d))

counts = Counter(d for *_, d in out)
with open("risk_mapping.tsv", "w") as fh:
    fh.write("run\texperiment_srx\tgsm\tdiagnosis\n")
    for run, exp, g, d in sorted(out):
        fh.write(f"{run}\t{exp}\t{g}\t{d}\n")
print(f"study: {study}; matrix samples: {len(gsm)}; runs mapped: {len(out)}")
print(f"diagnosis counts: {dict(counts)}")
