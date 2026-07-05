mkdir -p pdbs
while read id chain; do
  pdb_fetch "$id" \
    | pdb_selchain -"$chain" \
    | pdb_delhetatm \
    | pdb_tidy \
    > "pdbs/${id}_${chain}.pdb"
done < representatives.txt