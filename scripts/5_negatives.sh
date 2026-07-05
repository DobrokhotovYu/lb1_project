#   (reviewed:true) AND (fragment:false) NOT (xref:pfam-PF00014)

curl -s "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29+NOT+%28xref%3Apfam-PF00014%29" \
  > negatives_all.fasta
grep -c '>' negatives_all.fasta