#!/usr/bin/env python3
"""
plot_errors.py — Type I (false positive) and Type II (false negative) errors
as a function of the E-value threshold.

Reads the same files as the scorer:
  setN.ids : ground truth  (seq_id <TAB> 1|0)
  setN.tbl : hmmsearch --tblout output (hits only)

Counts are summed across both CV folds, so each curve is the total error count
over the entire evaluation set at each threshold. Marks the threshold that
minimises total errors (FP + FN). Saves PNG + PDF.

"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

FOLDS = ["set1", "set2"]
THRESHOLDS = [10 ** (-e) for e in range(1, 13)]   # 1e-1 ... 1e-12


def load_labels(path):
    d = {}
    for line in open(path):
        sid, lab = line.rstrip("\n").split("\t")
        d[sid] = int(lab)
    return d


def load_evalues(path):
    ev = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.split()
        if len(c) < 5:
            continue
        sid, e = c[0], float(c[4])
        if sid not in ev or e < ev[sid]:
            ev[sid] = e
    return ev


def counts(labels, ev, thr):
    FP = FN = 0
    for sid, truth in labels.items():
        pred = 1 if (ev.get(sid) is not None and ev[sid] <= thr) else 0
        if   truth == 0 and pred == 1: FP += 1
        elif truth == 1 and pred == 0: FN += 1
    return FP, FN


def main():
    data = [(load_labels(f"{s}.ids"), load_evalues(f"{s}.tbl")) for s in FOLDS]

    fps, fns = [], []
    for thr in THRESHOLDS:
        fp = fn = 0
        for labels, ev in data:
            a, b = counts(labels, ev, thr)
            fp += a
            fn += b
        fps.append(fp)
        fns.append(fn)

    # threshold minimising total errors (FP + FN)
    totals = [fp + fn for fp, fn in zip(fps, fns)]
    best_i = min(range(len(totals)), key=lambda i: totals[i])
    best_thr = THRESHOLDS[best_i]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(THRESHOLDS, fps, "-o", color="#1f77b4", label="False positives (Type I)")
    ax.plot(THRESHOLDS, fns, "-o", color="#ff7f0e", label="False negatives (Type II)")
    ax.axvline(best_thr, color="#d62728", ls="--", lw=1.3,
               label=f"Min-error threshold = {best_thr:.0e} "
                     f"(FP {fps[best_i]}, FN {fns[best_i]})")

    ax.set_xscale("log")
    ax.set_xlabel("E-value threshold (log scale)")
    ax.set_ylabel("Count")
    ax.set_title("Type I and Type II errors as a function of the E-value threshold")
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=13))
    ax.set_ylim(bottom=-0.5)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="upper center", framealpha=0.95)
    fig.tight_layout()
    fig.savefig("errors_vs_threshold.png", dpi=200)
    fig.savefig("errors_vs_threshold.pdf")
    print(f"wrote errors_vs_threshold.png / .pdf  "
          f"(min total errors {totals[best_i]} at {best_thr:.0e})")


if __name__ == "__main__":
    main()
