use the mipgen pipeline builder to analyze a set or to form a template for multiple analyses

edit the pipeline builder to change your desired default behavior
alternatively, edit the output of the sh file to customize analysis
see the options for each of the mipgen scripts by providing the -h flag

mipgen_smmip_collapser outputs the performance of each sample, each mip across all samples, and each mip for each sample

other tools can be used to:
- facilitate assessing mip design quality by building a UCSC track (generate_ucsc_track.py)
- pull out exons from a refGene file from a list of gene symbols (extract_coding_gene_exons.sh)
- assess basewise coverage from mip analysis output (complexity_to_basewise.sh)
- print a fasta of mip targets for specialized alignment purposes (prepare_mip_reference.py)

the steps of analysis are generally as follows:

step 0:
convert qseq format to (gzipped) fastq format with mipgen_qseq_to_fastq

step 1:
demultiplex samples by either inserting index into header with mipgen_fq_cutter or splitting runs into one file per sample

-optionally select specific index sequences (with or without 1bp Hamming distance)

-for dual indexed reads, paste indices together like so:

```
  paste index1.fq index2.fq | awk '(NR-1)%2==0 {print $1} (NR-1)%2==1 {print $1 $2}' > new_index.fq
```

step 2:
use PEAR to merge paired MIP reads with overlap

http://bioinformatics.oxfordjournals.org/content/early/2013/11/10/bioinformatics.btt593.long

step 3:
finish preprocessing reads with mipgen_fq_cutter

-move smmip tag to the header

-optionally trim reads

step 4:
map with bwa (mem) to human reference: Shendure's at /net/shendure/vol10/nobackup/shared/alignments/bwa-0.6.1/human_g1k_hs37d5/hs37d5.fa

-ensure version of bwa is above 0.6.1

step 5:
collapse smmip TDRGs with mipgen_smmip_collapser

-specify output directory and prefix and tag length

