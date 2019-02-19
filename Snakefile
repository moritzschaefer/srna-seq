import pandas as pd
from snakemake.utils import validate, min_version

include: "../../lib/srna-seq-star-deseq2/Snakefile"

##### load config and sample sheets #####

configfile: "config.yaml"
validate(
    config, schema="../../lib/srna-seq-star-deseq2/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"]).set_index("sample", drop=False)
validate(
    samples,
    schema="../../lib/srna-seq-star-deseq2/schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(
    units, schema="../../lib/srna-seq-star-deseq2/schemas/units.schema.yaml")

##### target rules #####

rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg"

##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### setup report #####

report: "report/workflow.rst"

include: "rules/rna_type_comparison.smk"
include: "rules/diffexp.smk"
