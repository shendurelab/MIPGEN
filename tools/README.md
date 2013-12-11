step 1:
demultiplex samples by either inserting index into header with mipgen_fq_cutter or splitting runs into one file per sample
-optionally select specific index sequences (with or without 1bp Hamming distance)
-for dual indexed reads, paste indices together like so:

```
  paste index1.fq index2.fq | awk '(NR-1)%2==0 {print $1} (NR-1)%2==1 {print $1 $2}' > new_index.fq

step 2:
use PEAR to merge paired MIP reads with overlap
http://bioinformatics.oxfordjournals.org/content/early/2013/11/10/bioinformatics.btt593.long
step 3:
finish preprocessing reads with mipgen_fq_cutter
-move smmip tag to the header
-optionally trim reads
step 4:
map with bwa (mem) to human reference at /net/shendure/vol10/nobackup/shared/alignments/bwa-0.6.1/human_g1k_hs37d5/hs37d5.fa
-ensure version of bwa is above 0.6.1
step 5:
collapse smmip TDRGs with mipgen_smmip_collapser
-specify output directory and prefix and tag length

