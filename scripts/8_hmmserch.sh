#!/usr/bin/env bash
# 8_hmmsearch.sh — searches the HMM against each fold and build ground-truth labels
set -euo pipefail
cd ~/Downloads/lb1_project

HMM=kunitz_kd.hmm

for f in 1 2; do
    pos="pos_folds/stdin.part_00${f}.fasta"
    neg="neg_folds/stdin.part_00${f}.fasta"

    # 1. assembles the evaluation set (positives + negatives for this fold)
    cat "$pos" "$neg" > "set${f}.fasta"

    # 2. searches the model (--max = full sensitivity; --tblout = per-hit table)
    hmmsearch --max --tblout "set${f}.tbl" "$HMM" "set${f}.fasta" > "set${f}.log"

    # 3. ground-truth labels: every sequence ID with its class (1 = Kunitz, 0 = not)
    grep '^>' "$pos" | sed 's/^>//; s/ .*//' | awk '{print $1"\t1"}'  > "set${f}.ids"
    grep '^>' "$neg" | sed 's/^>//; s/ .*//' | awk '{print $1"\t0"}' >> "set${f}.ids"

    # 4. report
    searched=$(grep -c '>' "set${f}.fasta")
    hits=$(grep -vc '^#' "set${f}.tbl")
    labeled=$(wc -l < "set${f}.ids")
    npos=$(awk '$2==1' "set${f}.ids" | wc -l)
    echo "set${f}: ${searched} searched | ${hits} in table | ${labeled} labeled (${npos} positive)"
done

# 5. format check: table IDs must match label IDs (both should be sp|ACC|NAME)
echo "--- ID format check (should match) ---"
echo "labels:"; head -3 set1.ids
echo "table :"; grep -v '^#' set1.tbl | awk '{print $1}' | head -3