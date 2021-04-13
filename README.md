# Usage

This pipeline depends on a submodule, so you need to run
git submodule init && git submodule update

sRNA-seq pipeline for read count- and DEG-analysis

There are three analyses:
- "normal" miRNA diffexp of RNAi mutants (in order to identify MI miRNAs)
- tRNA fragments
- AGO RIP-seq and Ago1-overexpression (RIP-seq is in ../rip-seq/)


In order to add an automatic quality control step (multiqc) just add qc/multiqc.html to the pipeline

*special* normalization means normalized using tRNA,snRNA and snoRNAs
