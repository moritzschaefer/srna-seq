import pandas as pd
from snakemake.utils import validate, min_version

##### load config and sample sheets #####

validate(
    config, schema="schemas/config.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(
    units, schema="schemas/units.schema.yaml")

include: "./lib/srna-seq-star-deseq2/Snakefile"
##### target rules #####

# this should go into config.yaml
RNA_TYPES = ['tdrs', 'mirbase_mmu21']

rule all:
    input:
        expand(["results/diffexp.{rna_type}/{contrast}.diffexp.tsv",
                "results/diffexp.{rna_type}/{contrast}.ma-plot.svg"],
               rna_type=RNA_TYPES,
               contrast=config["diffexp"]["contrasts"]),
        expand("counts/{rna_type}.cpm.tsv", rna_type=RNA_TYPES),
        "counts/wen15_mirtrons.tsv",  # the unmapped (mirbase)  mirtrons from wen15
        "counts/wen15_mirtrons.cpm.tsv",
        expand("plot/{rna_type}.pca.svg", rna_type=RNA_TYPES),
        
        "plot/expression_comparison.png",  # compare trnas, snornas, srnas, mirnas (to see what fits for normalization)
        "qc/multiqc.html",  # QC of trimmed reads 

        # only for tRNAs:
        "plot/tdrs_pca.png",  # TODO delete? (because deseq plots it...)
        "plot/tdrs_correlation.png",

##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### setup report #####

report: "report/workflow.rst"

include: "rules/wen_annotation.smk"
include: "rules/rna_type_comparison.smk"
include: "rules/diffexp.smk"
include: "rules/trnas.smk"
