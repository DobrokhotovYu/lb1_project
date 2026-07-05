curl -s "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28xref%3Apfam-PF00014%29+AND+%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29" \
  > positives_all.fasta
grep -c '>' positives_all.fasta