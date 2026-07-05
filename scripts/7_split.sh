# shuffle + halve the positives (353 -> ~176 / ~177)
seqkit shuffle -s 42 positives.fasta | seqkit split2 -p 2 -O pos_folds -f
# shuffle + halve the negatives
seqkit shuffle -s 42 negatives_all.fasta | seqkit split2 -p 2 -O neg_folds -f

ls pos_folds neg_folds
grep -c '>' pos_folds/*.fasta neg_folds/*.fasta