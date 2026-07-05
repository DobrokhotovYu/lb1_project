# HMM-Based Detection of BPTI/Kunitz Domains from Structure-Derived Alignments

![University of Bologna](https://img.shields.io/badge/university-Bologna-red.svg)
![Bioinformatics Lab 1](https://img.shields.io/badge/course-Bioinformatics%20Lab%201-blueviolet)
![HMMER](https://img.shields.io/badge/tool-HMMER-yellow)
![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue)
![CD-HIT](https://img.shields.io/badge/tool-CD--HIT-orange)
![MMseqs2](https://img.shields.io/badge/tool-MMseqs2-green)
![PDBeFold](https://img.shields.io/badge/tool-PDBeFold-lightblue)
![Project Status](https://img.shields.io/badge/status-complete-success)

This repository contains the complete pipeline and materials for building, evaluating, and validating a profile HMM that detects **BPTI/Kunitz-type protease-inhibitor domains** in protein sequences, produced for the final assessment of **Laboratory of Bioinformatics 1 @ University of Bologna**.

---

## Project Overview

The Kunitz domain is a compact (~58-residue), cysteine-rich serine-protease-inhibitor motif stabilised by three conserved disulfide bonds. This project builds a **structure-informed profile HMM** from a multiple structural alignment of curated, high-resolution Kunitz domains, then evaluates it against UniProtKB/Swiss-Prot with 2-fold cross-validation and scans the whole database to characterise the domain's distribution.

**Headline results:** cross-validated **MCC = 0.9929** (mean of the two folds); a whole-Swiss-Prot scan recovered **364 / 368** annotated Kunitz proteins with **zero false positives** at an E-value threshold of 1e-3.

---

## 🗂️ Repository Structure

```
lb1_project/
├── scripts/                                    # pipeline scripts
│   ├── 1_import_csv.py                          #   RCSB CSV report -> kunitz.fasta
│   ├── 2_select_representatives.py             #   length-gate + best-resolution selection
│   ├── 3_fetch_chains.sh                        #   extract single Kunitz chains
│   ├── 4_positives.sh   5_negatives.sh          #   UniProt positive / negative sets
│   ├── 6_leakage_removal.sh                     #   remove training proteins from positives
│   ├── 7_split.sh                               #   2-fold cross-validation split
│   ├── 8_hmmsearch.sh                           #   hmmsearch per fold + ground-truth labels
│   ├── 9_score.py                               #   confusion matrix, E-value sweep, CV MCC
│   ├── plot_mcc.py   plot_errors.py             #   result figures
│   ├── 10_swissprot_scan.sh                     #   hmmsearch vs all Swiss-Prot
│   ├── 11_swissprot_scan.py                     #   contingency table + distribution
│   └── 12_superkunitz.sh                        #   extreme multi-domain proteins
├── clustering/                                  # CD-HIT / MMseqs2 outputs
├── pdbs/                                        # extracted single-chain structures
├── pos_folds/   neg_folds/                      # 2-fold split FASTA parts
├── tmp/                                         # scratch (MMseqs2 temp)
│
├── rcsb_pdb_custom_report_20260704202953.csv    # raw RCSB report
├── kunitz.fasta                                 # all 135 candidate chains
├── representatives.{txt,fasta}                  # 12 selected representatives
├── representatives.log.tsv                      # selection decision log
├── PDBeFold_alignment.fasta                     # raw PDBeFold multiple alignment
├── seed.fasta                                   # curated 9-domain seed alignment
├── seed_seqs.fasta                              # ungapped seed sequences (for BLAST)
├── kunitz_kd.hmm                                # trained profile HMM
├── selfcheck.out                                # consistency-test output
│
├── positives_all.fasta                          # 368 PF00014 proteins
├── positives.fasta                              # 353 after leakage removal
├── negatives_all.fasta.part-a{a,b}              # negatives (split for GitHub size limit)
├── seed_vs_pos.tsv   training_ids.txt           # leakage-removal intermediates
├── posdb.*                                      # BLAST database (leakage removal)
│
├── set1.{ids,tbl,log,fasta.gz}                  # CV fold 1
├── set2.{ids,tbl,log,fasta.gz}                  # CV fold 2
│
├── swissprot_all.fasta.part-a{a,b}              # whole-DB scan input (split for GitHub size limit)
├── swissprot_scan.{tbl,log}                     # whole-DB scan output
├── unannotated_hits.tsv                         # predicted-but-unannotated hits
│
├── mcc_vs_threshold.{pdf,png}                   # figures
├── errors_vs_threshold.{pdf,png}
├── kunitz_score_distribution.{pdf,png}
│
├── environment.yml                              # conda environment
├── LICENSE                                      # license
└── README.md
```

---

## Project Workflow Summary

### 1. Data Acquisition and Preprocessing

Structures were retrieved from **RCSB PDB** with the query:

- Pfam ID: `PF00014`
- Resolution ≤ 3 Å
- Sequence length 45 ≤ length ≤ 80

A custom report was exported as CSV (Entry ID, Entity ID, Auth Asym ID, Sequence, Refinement Resolution, Annotation Identifier). Sequences were extracted with `scripts/1_import_csv.py`, yielding **135 chains**.

### 2. Redundancy Reduction

Clustered with **CD-HIT v4.8.1** (primary) and cross-checked with **MMseqs2**:

```bash
cd-hit -i kunitz.fasta -o kunitz_nr.fasta -c 0.90 -n 5 -d 0 -M 0 -T 0
mmseqs easy-cluster kunitz.fasta clusterRes tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
```

> Redundancy reduction was performed independently with CD-HIT (90% identity) and MMseqs2 (90% identity, 80% coverage). After selecting the highest-resolution, canonical-length domain per cluster, both approaches yielded concordant non-redundant seed sets (11/12 representatives identical), indicating the seed is robust to the choice of clustering method.

### 3. Representative Selection

CD-HIT's default (longest-sequence) representative is inadequate for a structure-based model, as it favours length outliers over clean domains. Representatives were reselected per cluster by **canonical length gate (52–66 aa) then best resolution** using `scripts/2_select_representatives.py`, then hand-reviewed.

### 4. Chain Extraction

Representative chains were fetched with `scripts/3_fetch_chains.sh`, followed by a length sanity check and manual review to confirm the correct chain was selected:

```bash
for f in pdbs/*.pdb; do
  echo "$f -> $(grep -c '^ATOM.* CA ' "$f") CA atoms"
done
```

### 5. Structural Alignment

A multiple structural alignment of the representative domains was produced with **PDBeFold**:

| Metric | Value |
| --- | --- |
| Aligned residues | 52 |
| Overall RMSD | 0.878 Å |
| Aligned SSEs | 3 |
| Overall Q-score | 0.732 |

> Inspection of the structural alignment revealed three engineered reduced-disulfide BPTI variants (`3wny`, `5jb7`, `1yld`) carrying Cys→Ala/Abu substitutions at conserved positions; these were excluded from the seed, yielding **9 native domains with all six catalytic cysteines fully conserved**.

### 6. HMM Construction

```bash
hmmbuild --amino --informat afa kunitz_kd.hmm seed.fasta
```

Consistency test — the model recognises its own seed:

```bash
hmmsearch --max kunitz_kd.hmm representatives.fasta > selfcheck.out
```

### 7. Validation Datasets

Positive and negative sets were retrieved from UniProtKB/Swiss-Prot (`scripts/4_positives.sh`, `scripts/5_negatives.sh`):

- **Positives:** reviewed, `PF00014`, non-fragment (368 proteins)
- **Negatives:** reviewed, non-fragment, non-`PF00014` (~566k proteins)

Training proteins and close homologs were removed from the positive set with `scripts/6_leakage_removal.sh` (BLAST, ≥95% identity, ≥50 aligned residues) — **15 entries removed**, leaving 353.

### 8. Cross-Validation and Scoring

Both classes were split into two folds (`scripts/7_split.sh`, fixed seed), searched with the HMM (`scripts/8_hmmsearch.sh`), and scored across E-value thresholds 1e-1 … 1e-12 with `scripts/9_score.py`. The MCC-optimal threshold chosen on one fold was evaluated on the other (and vice versa). Figures were produced with `plot_mcc.py` and `plot_errors.py`.

### 9. Whole-Swiss-Prot Scan

The final model was searched against all reviewed Swiss-Prot (`scripts/10_swissprot_scan.sh`) and analysed (`scripts/11_swissprot_scan.py`, `scripts/12_superkunitz.sh`) to build the predicted-vs-annotated contingency table and characterise the domain's taxonomic and multi-domain distribution.

---

## 🛠️ Environment Setup

```bash
conda env create -f environment.yml
conda activate kunitz
```

**Key dependencies:** Python 3.11+ (Matplotlib), HMMER, BLAST+, CD-HIT, MMseqs2, SeqKit, pdb-tools. Structural alignment used the PDBeFold web server.

---

## 📈 Performance Results

**2-fold cross-validation** (threshold chosen on the other fold):

| Fold evaluated | Threshold | TP | FP | FN | Sensitivity | Precision | MCC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| set2 (from set1) | 1e-2 | 175 | 1 | 1 | 0.9943 | 0.9943 | 0.9943 |
| set1 (from set2) | 1e-3 | 174 | 0 | 3 | 0.9831 | 1.0000 | 0.9915 |

> **Mean cross-validated MCC = 0.9929.** MCC is maximal and threshold-insensitive across a stable plateau (1e-2 – 1e-5).

**Whole-Swiss-Prot scan** (E ≤ 1e-3, 566,200 proteins): 364 / 368 annotated Kunitz proteins recovered, **0 false positives**.

**Missed known Kunitz domains** (false negatives — divergent invertebrate/nematode domains under-represented in the seed):

| UniProt ID | Entry | Note |
| --- | --- | --- |
| A0A1Q1NL17 | HA11_HYAAI | divergent (tick) Kunitz |
| Q8WPG5 | KUNI_ORNKA | divergent (tick) Kunitz |
| O62247 | BLI5_CAEEL | nematode BLI-5 |
| D3GGZ8 | BLI5_HAECO | nematode BLI-5 |

The Kunitz domain is distributed across a single-to-many-domain architecture (299 single-domain proteins up to papilins with 15–16 domains) and concentrated taxonomically in **venomous metazoa** (spiders, snakes, sea anemones) alongside **mammalian coagulation/protease-inhibitor proteins**.

---

## Author

**Iurii Dobrokhotov** 
## License

See [`LICENSE`](LICENSE). Data derived from RCSB PDB and UniProtKB/Swiss-Prot remains subject to their respective terms.
