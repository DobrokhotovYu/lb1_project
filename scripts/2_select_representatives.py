#!/usr/bin/env python3
"""
select_representatives.py

Pick one structural representative per CD-HIT cluster for a structure-informed
profile-HMM seed of the BPTI/Kunitz domain.

Rule (per cluster):
    1. LENGTH GATE : keep only members MIN_LEN..MAX_LEN residues (canonical ~58).
                     This drops "domain + extra chain" constructs that CD-HIT
                     otherwise promotes via its longest-sequence default.
    2. RESOLUTION  : among gated members, take the lowest (best) resolution.
    3. TIE-BREAK   : if resolutions tie, prefer length closest to CENTER.
Clusters with no canonical-length member are written to the log as FLAGGED
"""
import csv
import re

# ------------------------------------------------------------------
CLSTR      = "clustering/cd-hit clustering/kunitz_nr.fasta.clstr"                    # CD-HIT .clstr output
REPORT     = "rcsb_pdb_custom_report_20260704202953.csv"              # RCSB custom-report CSV
OUT_PREFIX = "representatives"                          # output file prefix
MIN_LEN, MAX_LEN, CENTER = 52, 66, 58                   # length gate + tie-break
# -----------------------------------------------------------------------------


def parse_resolution(cell):
    """RCSB may list several values ('2.15, 1.6'); takes the best (lowest)."""
    vals = []
    for part in cell.replace(";", ",").split(","):
        try:
            vals.append(float(part.strip()))
        except ValueError:
            pass
    return min(vals) if vals else None


def load_report(path):
    """RCSB custom-report CSV -> {id: (resolution, length, sequence)}.
    Handles the two-row grouped header and continuation rows (extra chains,
    empty first column). Columns located by header name, so order is flexible.
    """
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))

    header_idx = None
    for i, row in enumerate(rows[:5]):
        joined = " | ".join(c.lower() for c in row)
        if "sequence" in joined and "auth asym" in joined:
            header_idx = i
            break
    if header_idx is None:
        raise SystemExit("ERROR: header row with 'Sequence' + 'Auth Asym ID' "
                         "not found in " + path)

    header = [c.strip().lower() for c in rows[header_idx]]

    def col(*needles):
        for j, name in enumerate(header):
            if all(n in name for n in needles):
                return j
        raise SystemExit(f"ERROR: column matching {needles} not found.")

    c_entity, c_res, c_seq, c_chain = (col("entity id"), col("resolution"),
                                       col("sequence"), col("auth asym"))
    meta = {}
    for row in rows[header_idx + 1:]:
        if not row or not row[c_entity].strip():
            continue
        entity = row[c_entity].strip()
        res    = parse_resolution(row[c_res])
        seq    = row[c_seq].strip().replace(" ", "")
        chain  = row[c_chain].strip().split(",")[0].strip()
        if not seq or res is None:
            continue
        rid = f"{entity.split('_')[0].lower()}_{chain}"
        meta.setdefault(rid, (res, len(seq), seq))
    return meta


def load_clusters(path):
    """CD-HIT .clstr -> list of member-id lists."""
    clusters, cur = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                cur = []
                clusters.append(cur)
            else:
                m = re.search(r">(\S+?)\.\.\.", line)
                if m and cur is not None:
                    cur.append(m.group(1))
    return clusters


def main():
    meta = load_report(REPORT)
    clusters = load_clusters(CLSTR)

    missing = {i for cl in clusters for i in cl} - set(meta)
    if missing:
        print(f"WARNING: {len(missing)} clustered id(s) not in report "
              f"(id-format mismatch?): {sorted(missing)[:5]}")

    clean, flagged, log = [], [], []
    for k, cl in enumerate(clusters):
        ids = [i for i in cl if i in meta]
        if not ids:
            continue
        longest = max(ids, key=lambda i: (meta[i][1], -meta[i][0]))
        canon = [i for i in ids if MIN_LEN <= meta[i][1] <= MAX_LEN]
        if canon:
            rep = min(canon, key=lambda i: (meta[i][0], abs(meta[i][1] - CENTER)))
            clean.append(rep)
            status = "clean"
        else:
            rep = min(ids, key=lambda i: meta[i][0])
            flagged.append(rep)
            status = "FLAGGED_no_canonical_member"
        res, length, _ = meta[rep]
        log.append(dict(cluster=k, size=len(ids), representative=rep,
                        resolution=res, length=length, default_longest=longest,
                        override=("yes" if longest != rep else "no"),
                        status=status))

    selected = sorted(clean, key=lambda i: meta[i][0])

    with open(f"{OUT_PREFIX}.txt", "w") as f:
        for i in selected:
            pdb, chain = i.split("_")
            f.write(f"{pdb} {chain}\n")
    with open(f"{OUT_PREFIX}.fasta", "w") as f:
        for i in selected:
            f.write(f">{i}\n{meta[i][2]}\n")
    with open(f"{OUT_PREFIX}.log.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(log[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(log)

    print(f"clusters processed   : {len(log)}")
    print(f"clean representatives : {len(clean)}  -> {OUT_PREFIX}.txt / .fasta")
    print(f"flagged (manual review): {len(flagged)}  {flagged if flagged else ''}")


if __name__ == "__main__":
    main()
