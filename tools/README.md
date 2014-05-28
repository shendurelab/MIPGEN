Use the mipgen pipeline builder to analyze a set or to form a template for multiple analyses

```
python mipgen_pipeline_builder.py
```

Edit the pipeline builder to change your desired default behavior

Alternatively, edit the output of the sh file to customize analysis. Analysis shell scripts link together the other scripts in the tools directory along with alignment with BWA to lead to a SAM file ready for variant calling. Memory text files record your answers to the builder script to clarify how the analysis script was composed. Feeding memory files into the builder script (with or without modification) enables fast rebuilding of the analysis script.

See the options for each of the mipgen scripts by providing the -h flag

mipgen_smmip_collapser outputs the performance of each sample, each mip across all samples, and each mip for each sample

Other tools can be used to:
- facilitate assessing mip design quality by building a UCSC track (generate_ucsc_track.py)
- pull out exons from a refGene file from a list of gene symbols (extract_coding_gene_exons.sh)
- assess basewise coverage from mip analysis output (complexity_to_basewise.sh)
- print a fasta of mip targets for specialized alignment purposes (prepare_mip_reference.py)

It is first necessary to compile the Cython module genome_sam_collapser.pyx into a .so file

```
python setup.py build_ext --inplace
```

The steps of analysis are generally as follows:

Step 0 (starting from qseq format ensures fastq headers are properly formatted):
convert qseq format to (gzipped) fastq format with mipgen_qseq_to_fastq

Step 1:
demultiplex samples by either inserting index into header with mipgen_fq_cutter or splitting runs into one file per sample

-optionally select specific index sequences (with or without 1bp Hamming distance)

-for dual indexed reads, paste indices together like so:

```
  paste index1.fq index2.fq | awk '(NR-1)%2==0 {print $1} (NR-1)%2==1 {print $1 $2}' > new_index.fq
```

Step 2:
use PEAR to merge paired MIP reads with overlap

http://bioinformatics.oxfordjournals.org/content/early/2013/11/10/bioinformatics.btt593.long

Step 3:
finish preprocessing reads with mipgen_fq_cutter

-move smmip tag to the header

-optionally trim reads

Step 4:
map with bwa (mem) to human reference

-ensure version of bwa is above 0.6.1

Step 5:
collapse smmip TDRGs with mipgen_smmip_collapser

-specify output directory and prefix and tag length

-disable collapsing (and partition reads) with -w option

Step 6:
review output files:

on target, deduped/smc-reads
-all_reads.unique.sam
-[sample].unique.sam

reads with softclipped MIP arms
-softclipped.sam

paired end reads mapping kilobases/chromosomes apart
-discordant_arms.sam

paired end reads missing their pairs
-unpaired_reads.sam

valid reads capturing unintended genomic loci
-off_target_reads.sam

reads flagged as improper alignments by BWA
-improper_pairs.sam

reads with unexpected/incorrectly formatted SAM fields
-strange_alignments.sam

reads with targeting arm sequencing not matching synthesized oligos
-imperfect_arms.sam

MIP summary information
-smMIPS_noPEAR.indexed.sort.collapse.complexity.txt
-smMIPS_noPEAR.indexed.sort.collapse.mipwise_summary.txt
-smMIPS_noPEAR.indexed.sort.collapse.samplewise_summary.txt
-smMIPS_noPEAR.indexed.sort.collapse.notes.txt

© University of Washington 2014
