# ungapped seed sequences (strips alignment gaps from seed.fasta)
seqkit seq -g seed.fasta > seed_seqs.fasta     # or: awk '/^>/{print;next}{gsub(/-/,"");print}'

makeblastdb -in positives_all.fasta -dbtype prot -out posdb
blastp -query seed_seqs.fasta -db posdb \
       -outfmt "6 sseqid pident length" -evalue 1e-3 > seed_vs_pos.tsv

# training-related hits: >=95% identity over >=50 aligned residues
awk '$2>=95 && $3>=50 {print $1}' seed_vs_pos.tsv | sort -u > training_ids.txt
wc -l training_ids.txt                          # the proteins that will be removed

# remove them from the positive set
seqkit grep -v -f training_ids.txt positives_all.fasta > positives.fasta
grep -c '>' positives.fasta                      #leakage-free positives