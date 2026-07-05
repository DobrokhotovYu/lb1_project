#!/usr/bin/env python3
"""
MCC as a function of the E-value threshold (per CV fold).

Reads the same files as the scoring script:
  setN.ids : ground truth  (seq_id <TAB> 1|0)
  setN.tbl : hmmsearch --tblout output (hits only)

Computes MCC across a range of E-value thresholds for each fold, plots them
on a log x-axis, and marks each fold's MCC-optimal threshold. Saves PNG + PDF.

"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

FOLDS = ["set1", "set2"]
THRESHOLDS = [10 ** (-e) for e in range(1, 13)]   # 1e-1 ... 1e-12
COLORS = {"set1": "#1f77b4", "set2": "#d62728"} 
PLATEAU = (1e-5, 1e-2)                             # stable operating region to shade


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


def mcc_at(labels, ev, thr):
    TP = FP = FN = TN = 0
    for sid, truth in labels.items():
        e = ev.get(sid)
        pred = 1 if (e is not None and e <= thr) else 0
        if   truth == 1 and pred == 1: TP += 1
        elif truth == 0 and pred == 1: FP += 1
        elif truth == 1 and pred == 0: FN += 1
        else:                          TN += 1
    denom = math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    return ((TP * TN - FP * FN) / denom) if denom else 0.0


def main():
    fig, ax = plt.subplots(figsize=(9, 5.5))
    # shaded stable operating region (plateau)
    ax.axvspan(PLATEAU[0], PLATEAU[1], color="#f0e68c", alpha=0.35, zorder=0,
               label=f"stable plateau ({PLATEAU[1]:.0e}\u2013{PLATEAU[0]:.0e})")
    all_mcc = []
    for s in FOLDS:
        labels, ev = load_labels(f"{s}.ids"), load_evalues(f"{s}.tbl")
        mccs = [mcc_at(labels, ev, t) for t in THRESHOLDS]
        all_mcc += mccs
        ax.plot(THRESHOLDS, mccs, "-o", color=COLORS.get(s, None),
                markersize=5, label=f"{s} (MCC)")
        # optimal threshold for this fold
        best_i = max(range(len(mccs)), key=lambda i: mccs[i])
        best_thr, best_mcc = THRESHOLDS[best_i], mccs[best_i]
        ax.axvline(best_thr, color=COLORS.get(s, "grey"), ls="--", lw=1.2, alpha=0.8)
        ax.annotate(f"{s} optimum\n{best_thr:.0e} (MCC {best_mcc:.4f})",
                    xy=(best_thr, best_mcc),
                    xytext=(6, -14 if s == "set2" else 6),
                    textcoords="offset points", fontsize=8,
                    color=COLORS.get(s, "black"))

    ax.set_xscale("log")
    ax.set_xlabel("E-value threshold (log scale)")
    ax.set_ylabel("Matthews Correlation Coefficient (MCC)")
    ax.set_title("MCC as a function of the E-value threshold")
    lo = min(all_mcc)
    ax.set_ylim(lo - 0.01, 1.001)
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=13))
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="lower center", framealpha=0.95)
    fig.tight_layout()
    fig.savefig("mcc_vs_threshold.png", dpi=200)
    fig.savefig("mcc_vs_threshold.pdf")
    print("wrote mcc_vs_threshold.png and mcc_vs_threshold.pdf")


if __name__ == "__main__":
    main()
