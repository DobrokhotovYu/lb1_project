#!/usr/bin/env python3
"""
9_score.py — evaluate the Kunitz HMM with 2-fold cross-validation.

For each fold it reads:
  setN.ids  : ground truth, one "seq_id<TAB>label" per line (1=Kunitz, 0=not)
  setN.tbl  : hmmsearch --tblout output (only sequences that HIT the model)

Key correctness point: hmmsearch lists only hits. Every sequence in setN.ids
that is NOT in setN.tbl (or is above the E-value threshold) is a predicted
NEGATIVE. That is how the ~283k true negatives are counted.

It sweeps the E-value threshold, builds a confusion matrix at each, computes
MCC / sensitivity / precision / accuracy, then does proper cross-validation:
pick the MCC-optimal threshold on one fold, report the OTHER fold at that
threshold (and vice versa) — so the operating point is never chosen on the
data it is scored on.

"""
import math
import re

FOLDS = ["set1", "set2"]
# E-value thresholds to scan
THRESHOLDS = [10 ** (-e) for e in range(1, 13)]


def load_labels(path):
    """seq_id -> 0/1 ground-truth label."""
    labels = {}
    with open(path) as fh:
        for line in fh:
            sid, lab = line.rstrip("\n").split("\t")
            labels[sid] = int(lab)
    return labels


def load_evalues(path):
    """seq_id -> best (smallest) full-sequence E-value from hmmsearch --tblout.

    tblout columns (whitespace-separated):
      1 target_name  2 target_acc  3 query_name  4 query_acc  5 full_E-value ...
    """
    ev = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            sid = cols[0]
            evalue = float(cols[4])
            if sid not in ev or evalue < ev[sid]:
                ev[sid] = evalue
    return ev


def confusion(labels, evalues, thr):
    """Predicted positive if sequence is in the table AND E-value <= thr."""
    TP = FP = FN = TN = 0
    for sid, truth in labels.items():
        e = evalues.get(sid)
        pred = 1 if (e is not None and e <= thr) else 0
        if   truth == 1 and pred == 1: TP += 1
        elif truth == 0 and pred == 1: FP += 1
        elif truth == 1 and pred == 0: FN += 1
        else:                          TN += 1
    return TP, FP, FN, TN


def metrics(TP, FP, FN, TN):
    n = TP + TN + FP + FN
    acc = (TP + TN) / n if n else 0.0
    sen = TP / (TP + FN) if (TP + FN) else 0.0          # recall / TPR
    pre = TP / (TP + FP) if (TP + FP) else 0.0          # PPV
    denom = math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc = ((TP * TN - FP * FN) / denom) if denom else 0.0
    return acc, sen, pre, mcc


def main():
    data = {}
    for s in FOLDS:
        data[s] = (load_labels(f"{s}.ids"), load_evalues(f"{s}.tbl"))

    # --- full sweep per fold ---
    best = {}
    for s in FOLDS:
        labels, ev = data[s]
        print(f"\n=== {s}: E-value threshold sweep "
              f"({sum(labels.values())} pos / {len(labels)} total) ===")
        print(f"{'threshold':>10} {'TP':>5} {'FP':>5} {'FN':>5} {'TN':>8} "
              f"{'sens':>6} {'prec':>6} {'MCC':>7}")
        best_mcc, best_thr, best_cm = -2, None, None
        for thr in THRESHOLDS:
            cm = confusion(labels, ev, thr)
            acc, sen, pre, mcc = metrics(*cm)
            TP, FP, FN, TN = cm
            print(f"{thr:>10.0e} {TP:>5} {FP:>5} {FN:>5} {TN:>8} "
                  f"{sen:>6.3f} {pre:>6.3f} {mcc:>7.4f}")
            if mcc > best_mcc:
                best_mcc, best_thr, best_cm = mcc, thr, cm
        best[s] = (best_thr, best_mcc, best_cm)
        print(f"  -> best on {s}: threshold {best_thr:.0e}, MCC {best_mcc:.4f}")

    # --- 2-fold cross-validation: choose thr on one fold, report the other ---
    print("\n" + "=" * 60)
    print("2-FOLD CROSS-VALIDATION (threshold chosen on the OTHER fold)")
    print("=" * 60)
    pairs = [("set1", "set2"), ("set2", "set1")]
    cv_mccs = []
    for train, test in pairs:
        thr = best[train][0]                     # threshold optimised on `train`
        labels, ev = data[test]
        cm = confusion(labels, ev, thr)
        acc, sen, pre, mcc = metrics(*cm)
        TP, FP, FN, TN = cm
        cv_mccs.append(mcc)
        print(f"\nthreshold {thr:.0e} (from {train})  ->  evaluate on {test}:")
        print(f"  TP={TP}  FP={FP}  FN={FN}  TN={TN}")
        print(f"  sensitivity={sen:.4f}  precision={pre:.4f}  "
              f"accuracy={acc:.6f}  MCC={mcc:.4f}")
        if FP:
            fps = [sid for sid, t in labels.items()
                   if t == 0 and (ev.get(sid) is not None and ev[sid] <= thr)]
            print(f"  false positives ({len(fps)}): {', '.join(fps[:10])}"
                  + (" ..." if len(fps) > 10 else ""))
        if FN:
            fns = [sid for sid, t in labels.items()
                   if t == 1 and not (ev.get(sid) is not None and ev[sid] <= thr)]
            print(f"  false negatives ({len(fns)}): {', '.join(fns[:10])}"
                  + (" ..." if len(fns) > 10 else ""))

    print(f"\nCross-validated MCC: {cv_mccs[0]:.4f}, {cv_mccs[1]:.4f}  "
          f"(mean {sum(cv_mccs)/2:.4f})")


if __name__ == "__main__":
    main()
