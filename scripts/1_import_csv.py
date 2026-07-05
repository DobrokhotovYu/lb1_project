import csv

seen = set()
infile  = "rcsb_pdb_custom_report_20260704202953.csv"
outfile = "kunitz.fasta"

with open(infile) as f, open(outfile, "w") as out:
    rows = list(csv.reader(f))
    for line in rows[2:]:                       # skip the 2 grouped-header rows
        if not line or not line[0].strip():     # skip continuation rows (extra chains)
            continue
        ent   = line[0].strip()                 # e.g. 1AAL_1
        seq   = line[2].strip().replace(" ", "")
        chain = line[3].strip().split(",")[0].strip()  # first chain only
        annot = line[4] if len(line) > 4 else ""
        if not seq or "PF00014" not in annot:   # QC: keep only real PF00014 hits
            continue
        pdb = ent.split("_")[0]
        key = f"{pdb.lower()}_{chain}"          # header like >1aal_A
        if key in seen:                         # dedup identical (pdb, chain)
            continue
        seen.add(key)
        out.write(f">{key}\n{seq}\n")