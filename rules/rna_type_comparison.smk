# Here we compare and plot the relative expression of different RNA types which helps to select appropriate normalization factors

rule compare_rna_type_expression:
    input:
        mirbase_counts="counts/mirbase_mmu21.tsv",
        trna_counts="counts/mm10-tRNAs.tsv",
        gencode_counts="counts/gencode.vM20.annotation.tsv",
        gencode_annotation="ref/gencode.vM20.annotation.gff3"
    output:
        plot_path="plot/expression_comparison.png"
    conda:
        "../envs/pandas.yaml"
    params:
        samples=units.index.get_level_values('sample').tolist()
    script:
        "../scripts/compare_rna_type_expression.py"
