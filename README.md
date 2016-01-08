# Description
According to the beginning and the end conservative gene of gene cluster, extract the complete gene cluster sequences.


# Usage
perl geneClusterExtract.pl query_seq db_seq out_prefix extend_length

db_seq format is:  
\>start_gene direction  
fasta_seq  
\>end_gene direction  
fasta_seq  

the start_gene seq must at the beginning of db_seq, and direction must be 1 or -1  
extend_length default 0



Copyright (c) Zhihao Xie, xiezhihao1122@outlook.com
