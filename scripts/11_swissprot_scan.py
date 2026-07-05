#!/usr/bin/env python3
"""
Inputs (in the current directory):
  swissprot_scan.tbl   : hmmsearch --tblout of the HMM vs ALL reviewed Swiss-Prot
  positives_all.fasta  : the PF00014-annotated proteins (defines "annotated +")
  negatives_all.fasta  : the non-PF00014 proteins (only used for the DB total)

Produces:
  - a predicted x annotated 2x2 contingency table (printed)
  - unannotated_hits.tsv : proteins the model predicts Kunitz but that lack the
    PF00014 annotation, with an optional UniProt annotation-flag lookup
  - a taxonomy summary and a domains-per-protein summary (printed)
  - kunitz_score_distribution.png/.pdf : bit-score histogram, annotated vs not

"""
import re
import urllib.request
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCAN        = "swissprot_scan.tbl"
POS_FASTA   = "positives_all.fasta"
NEG_FASTA   = "negatives_all.fasta"
THRESHOLD   = 1e-3            # operating threshold (from the plateau)
CHECK_UNIPROT = True          # can be set to False to skip network lookups
OUT_HITS    = "unannotated_hits.tsv"
FIG         = "kunitz_score_distribution"


def count_fasta(path):
    try:
        return sum(1 for l in open(path) if l.startswith(">"))
    except FileNotFoundError:
        return None


def annotated_ids(path):
    ids = set()
    for l in open(path):
        if l.startswith(">"):
            ids.add(l[1:].split()[0])     # sp|ACC|NAME
    return ids


def uniprot_flags(seqid):
    """Look up whether an unannotated hit has Kunitz-consistent UniProt evidence."""
    acc = seqid.split("|")[1] if "|" in seqid else seqid
    url = (f"https://rest.uniprot.org/uniprotkb/{acc}"
           f"?format=tsv&fields=keyword,ft_domain,protein_families,cc_function")
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            lines = r.read().decode().splitlines()
        if len(lines) < 2:
            return "no_data"
        row = lines[1].lower()
        ev = []
        if "kunitz" in row:                                    ev.append("Kunitz")
        if "bpti" in row:                                      ev.append("BPTI")
        if "protease inhibitor" in row or "proteinase inhibitor" in row:
            ev.append("protease-inhibitor")
        return ";".join(ev) if ev else "none"
    except Exception:
        return "lookup_failed"


def main():
    annotated = annotated_ids(POS_FASTA)
    n_pos = len(annotated)
    n_neg = count_fasta(NEG_FASTA)
    total_db = (n_pos + n_neg) if n_neg is not None else None

    # parse the scan
    hits = {}
    for line in open(SCAN):
        if line.startswith("#"):
            continue
        c = line.split()
        if len(c) < 16:
            continue
        sid = c[0]
        full_e, full_score = float(c[4]), float(c[5])
        ndom = int(c[15])
        desc = " ".join(c[18:]) if len(c) > 18 else ""
        m = re.search(r"OS=(.+?)\s+OX=(\d+)", desc)
        org = m.group(1) if m else "?"
        if sid not in hits or full_e < hits[sid]["e"]:
            hits[sid] = dict(e=full_e, score=full_score, ndom=ndom, org=org)

    pred = {sid for sid, h in hits.items() if h["e"] <= THRESHOLD}

    TP = len(pred & annotated)
    FP = len(pred - annotated)
    FN = len(annotated - pred)
    TN = (total_db - TP - FP - FN) if total_db is not None else None

    print(f"Whole-Swiss-Prot scan  (threshold E <= {THRESHOLD:.0e})")
    print(f"  DB size: {total_db}   annotated PF00014: {n_pos}   hits in table: {len(hits)}\n")
    print("  Predicted x Annotated contingency:")
    print(f"                     annotated +   annotated -")
    print(f"    predicted +   {TP:>10}   {FP:>10}")
    tn_str = f"{TN}" if TN is not None else "(DB total unknown)"
    print(f"    predicted -   {FN:>10}   {tn_str:>10}")
    print(f"\n  agreement (TP): {TP}   |  missed known (FN): {FN}   "
          f"|  flagged-but-unannotated (FP): {FP}")

    # domains-per-protein among predicted positives
    ndoms = Counter(hits[s]["ndom"] for s in pred)
    print("\n  Domains per predicted-Kunitz protein:")
    for k in sorted(ndoms):
        print(f"    {k} domain(s): {ndoms[k]} proteins")

    # taxonomy of predicted positives
    orgs = Counter(hits[s]["org"] for s in pred)
    print("\n  Top organisms among predicted-Kunitz proteins:")
    for org, n in orgs.most_common(15):
        print(f"    {n:>4}  {org}")

    # unannotated hits (predicted + but not PF00014) -> the interesting set
    unann = sorted(pred - annotated, key=lambda s: -hits[s]["score"])
    print(f"\n  Writing {len(unann)} unannotated hits to {OUT_HITS}"
          + (" (with UniProt lookup)" if CHECK_UNIPROT else ""))
    with open(OUT_HITS, "w") as f:
        f.write("seq_id\tbit_score\tE_value\tn_domains\torganism\tuniprot_kunitz_evidence\n")
        for s in unann:
            h = hits[s]
            flag = uniprot_flags(s) if CHECK_UNIPROT else "not_checked"
            f.write(f"{s}\t{h['score']:.1f}\t{h['e']:.1e}\t{h['ndom']}\t{h['org']}\t{flag}\n")

    # score-distribution figure (annotated vs unannotated hits)
    ann_scores = [hits[s]["score"] for s in hits if s in annotated]
    un_scores  = [hits[s]["score"] for s in hits if s not in annotated]
    fig, ax = plt.subplots(figsize=(8, 5))
    bins = 40
    ax.hist(ann_scores, bins=bins, color="#1f77b4", alpha=0.7,
            label=f"PF00014-annotated ({len(ann_scores)})")
    ax.hist(un_scores, bins=bins, color="#d62728", alpha=0.8,
            label=f"not annotated ({len(un_scores)})")
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlabel("HMMER full-sequence bit-score")
    ax.set_ylabel("number of proteins")
    ax.set_title("Kunitz HMM hit-score distribution across Swiss-Prot")
    ax.set_yscale("log")           # log y so the sparse unannotated hits are visible
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{FIG}.png", dpi=200)
    fig.savefig(f"{FIG}.pdf")
    print(f"  wrote {FIG}.png / .pdf")


if __name__ == "__main__":
    main()
