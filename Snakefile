import pandas as pd
from snakemake.utils import validate, min_version

configfile: "config.yaml"

##### load config and sample sheets #####

validate(
    config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], dtype=str).set_index("sample", drop=False)
validate(
    samples,
    schema="schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(
    units, schema="../../lib/srna-seq-star-deseq2/schemas/units.schema.yaml")

include: "../../lib/srna-seq-star-deseq2/Snakefile"
##### target rules #####

rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["diffexp"]["contrasts"]),
        "counts/mirbase_mmu21.cpm.tsv",
        "counts/wen15_mirtrons.tsv",  # the unmapped (mirbase)  mirtrons from wen15
        "counts/wen15_mirtrons.cpm.tsv",
        # "results/pca.svg",
        "plot/expression_comparison.png"

##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### setup report #####

report: "report/workflow.rst"

include: "rules/wen_annotation.smk"
include: "rules/rna_type_comparison.smk"
include: "rules/diffexp.smk"
