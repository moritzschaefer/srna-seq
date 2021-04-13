import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler


# from moritzsphd.plot
def correlation_plot(df, scaling='', ax=None, cmap='Spectral'):
    '''
    Requires a normalized library in the first place (i.e. each sample should be CPM (or so) normalized)


    scaling: - 'standard': use standard scaler to scale the features
             - 'log2': use log2 eof all values
             - 'log2-standard' First apply log2, then standardize
             - '': don't scale
             'standard' and 'log2' help to avoid PCA domination by large genes

    '''
    data = df
    if 'log2' in scaling:
        data = np.log2(data + 1)

    if 'standard' in scaling:
        scaler = StandardScaler()

        # df is transposed because normally genes (corresponding to features) are the rows (but should be columns)
        data = scaler.fit_transform(data.T).T

    corr = data.corr()

    # Generate a mask for the upper triangle
    # mask = np.triu(np.ones_like(corr, dtype=np.bool))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    ax = sns.heatmap(corr, cmap=cmap, square=True, linewidths=.5, center=0,
                     cbar_kws={"shrink": .5}, ax=ax)

    return ax

df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0)

fig, ax = plt.subplots(figsize=(15, 15))
correlation_plot(df, ax=ax)

plt.tight_layout()
fig.savefig(snakemake.output[0])
