# Here we compare and plot the relative expression of different RNA types which helps to select appropriate normalization factors

rule compare_rna_type_expression:
    input:
        rna_type_counts=expand("counts/{annotation}.tsv", annotation=['mirbase_mmu21', 'mm10-tRNAs', 'gencode.vM20.annotation'])
    output:
        plot_path="plot/expression_comparison.png"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/compare_rna_type_expression.py"