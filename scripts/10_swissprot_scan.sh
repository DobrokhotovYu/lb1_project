cat positives_all.fasta negatives_all.fasta > swissprot_all.fasta
hmmsearch --max --tblout swissprot_scan.tbl kunitz_kd.hmm swissprot_all.fasta > swissprot_scan.log