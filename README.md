MIPGEN
======

One stop MIP design and analysis

Use MIPgen to design custom mip panels for target enrichment of moderate to high complexity DNA targets ranging from 120 to 250bp in size

To compile MIPgen, you will need a C++ compiler, such as GCC (http://gcc.gnu.org/). That's it! For running design and analysis, see below for dependencies

-----
BUILD MIPGEN
-----

To build MIPgen from source, simply enter the MIPGEN directory and type 'make'

The 'mipgen' executable can be used to perform designs

-----
MIPGEN PARAMETERS
-----

MIPgen accepts many options intended to maximize the flexibility of your designs

The only required parameters are the minimum and maximum sizes for the captured sequence, a prefix for output files and a BED file consisting of the coordinates of the desired targets, and the bwa- and samtools-indexed reference genome from which to pull the sequences

Other options provide control over tiling, scoring and oligo design behavior

MIPgen offers three scoring methods based on two methods of scoring: SVR and logistic regression. Both methods offer similar performance, with slight improvements seen with the SVR scoring. The 'mixed' scoring option will perform designs primarily with the faster logistic regression score and finish with the more accurate SVR scoring scheme.

Run ./mipgen -doc to see the full documentation

-----
MIPGEN PITFALLS
-----

-It is not recommended to include more than 8bp total of degenerate smMIP tags. Longer tags reduce specificity and complicate MIP rebalancing, leading to less complete coverage of targeted sites.

-Priority scores (used by default) modulate tiling behavior and density of tiling of low scoring regions. To reduce density, lower priority scores. In particular, for interrogation of very short loci (<10bp), priority scores should be set to 0 to prevent splitting the region into two MIP targets.

-Regions that do not map uniquely to the provided reference are not tiled by default. These and other bases that cannot be tiled will be reported in the <project_name>.coverage_failed.bed file.

-SVR scores below 1.4 and logistic scores below 0.7 rarely perform adequately. It is recommended to remove these probes from designed assays.

-----
OTHER NOTES
-----

-the subset of Boost required for compilation is included in the directory

-the makefile should compile the source without any external dependencies, but the executable that is generated makes system calls to SAMtools, BWA, Tabix, and, optionally, Tandem Repeats Finder

-the tools subdirectory has scripts for parsing fastqs prior to mapping and collapsing tag-defined read groups (TDRGs) and trimming targeting arms after mapping (Cython required)

-summary information is also printed as part of the collapsing process and requires NumPy and SciPy

-if there is at least 10bp of overlap between forward and reverse reads, it is recommended to merge read pairs into fr-reads using PEAR ( http://bioinformatics.oxfordjournals.org/content/early/2013/11/10/bioinformatics.btt593.long ) prior to processing and to utilize single end options
