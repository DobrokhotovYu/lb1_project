# pull extreme multi-domain proteins
awk -F'\t' 'NR>1' swissprot_scan.tbl >/dev/null   # (they're in the .tbl; find by ndom col 16)
grep -v '^#' swissprot_scan.tbl | awk '$16>=6 {print $1, $16, $0}' | awk '{print $1, $2}'