def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(units) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        target_counts="counts/{rna_type}.tsv",
        trna_counts="counts/mm10-tRNAs.tsv",
        snrna_counts="counts/snrnas.tsv",
        snorna_counts="counts/snornas.tsv",
        mt_trna_counts="counts/mt_trnas.tsv"
    output:
        all="deseq2/{rna_type}.all.rds",
        normal_normalized_counts="deseq2/{rna_type}.normal_normalized_counts.tsv",
        special_normalized_counts="deseq2/{rna_type}.special_normalized_counts.tsv"
    params:
        units=config["units"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.{rna_type}.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "deseq2/{rna_type}.all.rds"
    output:
        report("plot/{rna_type}.pca.svg", "../report/{rna_type}.pca.rst")
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.{rna_type}.log"
    script:
        "../scripts/plot-pca.R"


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]


rule deseq2:
    input:
        "deseq2/{rna_type}.all.rds"
    output:
        table=report("results/diffexp.{rna_type}/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp.{rna_type}/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{rna_type}.{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

