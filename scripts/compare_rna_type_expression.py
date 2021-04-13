import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

dfs = []
gencode_counts_df = pd.read_csv(
    snakemake.input.gencode_counts,
    sep='\t',
    index_col=0,
    usecols=['Geneid'] + snakemake.params.samples)

annotation = pd.read_csv(
    snakemake.input.gencode_annotation,
    comment='#',
    sep='\t',
    names=[
        "seq_id", "source", "type", "start", "end", "score", "strand", "phase",
        "attributes"
    ])

annotation['gene_id'] = annotation['attributes'].apply(
    lambda v: re.search('gene_id=([^;]*)', v).groups()[0])
annotation['gene_type'] = annotation['attributes'].apply(
    lambda v: re.search('gene_type=([^;]*)', v).groups()[0])

gencode_counts_df = gencode_counts_df.merge(
    annotation[['gene_id',
                'gene_type']].groupby('gene_id').first().reset_index(),
    left_on='Geneid',
    right_on='gene_id')

data = {}
for rna_type in ['snRNA', 'snoRNA', 'Mt_tRNA']:
    rna_type_sum = gencode_counts_df[gencode_counts_df['gene_type'] ==
                                     rna_type][snakemake.params.samples].sum(
                                         axis=0)
    data[rna_type] = rna_type_sum

data['miRNA'] = pd.read_csv(
    snakemake.input.mirbase_counts,
    sep='\t',
    index_col=0,
    usecols=['Geneid'] + snakemake.params.samples).sum(axis=0)
data['tRNA'] = pd.read_csv(
    snakemake.input.trna_counts,
    sep='\t',
    index_col=0,
    usecols=['Geneid'] + snakemake.params.samples).sum(axis=0)
data['tDR'] = pd.read_csv(
    snakemake.input.tdr_counts,
    sep='\t',
    index_col=0).sum(axis=0)
df = pd.DataFrame(data).T
df.index.name = 'rna_type'
fig, ax = plt.subplots(figsize=(14, 10))
sns.barplot(
    data=df.reset_index().melt(
        id_vars=['rna_type'],
        value_vars=snakemake.params.samples,
        var_name='sample',
        value_name='counts'),
    x='sample',
    y='counts',
    hue='rna_type',
    ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=20, ha='right')
ax.set_yscale('log')
plt.tight_layout()
plt.savefig(snakemake.output.plot_path, dpi=120)
