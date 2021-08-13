# Usage

This pipeline depends on a submodule, so you need to run
git submodule init && git submodule update

sRNA-seq pipeline for read count- and DEG-analysis of miRNA expression

In order to add an automatic quality control step (multiqc) just add qc/multiqc.html to the pipeline

*special* normalization means normalized using tRNA,snRNA and snoRNAs
