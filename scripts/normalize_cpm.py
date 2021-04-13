import re

import pandas as pd

df = pd.concat([pd.read_csv(input, sep='\t') for input in snakemake.input])
df['sample'] = df['file'].apply(lambda s: re.match('[^-]*', s).group())  # cut the unit (for downstream-compatibility reasons)
columns_order = df['sample'].drop_duplicates()
# we don't care about duplicate tRNAs (therefore we cut the last '-\d')
df['Geneid'] = df.apply(lambda row: row['tRNA'][:row['tRNA'].rfind('-')] + '_' + row['species'],
                     axis=1)  # tDR-name = Geneid, for compatibility with downstream analysis

# remove duplicates:
# - locMajorityCovg is ignored because it varies (slightly) across duplicates (weird..)
# - most rows with the same 'sample' and 'Geneid' also have the same values, so mean is same as drop_duplicates
# - only 21 duplicates have distinct (yet very similar) counts and percentages. We just mean them.

df = df[['sample', 'Geneid', 'counts']].groupby(['sample', 'Geneid']).mean().reset_index()

# all reads
df_pivot = df.pivot(index='Geneid', columns='sample', values='counts').fillna(0)
df_pivot[columns_order].astype(int).to_csv(snakemake.output['all'], sep='\t')

# CPM
df['cpm_counts'] = df[['sample', 'counts']].groupby('sample')['counts'].apply(lambda counts: 1e6 * counts / counts.sum())
df_pivot = df.pivot(index='Geneid', columns='sample', values='cpm_counts').fillna(0)
df_pivot[columns_order].to_csv(snakemake.output['all_cpm'], sep='\t')
