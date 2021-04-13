import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scikitplot.decomposition import plot_pca_2d_projection
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# from moritzsphd.plot
def pcaplot(df, dimension=2, scaling=None, ax=None, cmap='Spectral', label_dots=True):
    '''
    Requires a normalized library in the first place (i.e. each sample should be CPM (or so) normalized)

    df is transposed because normally genes (corresponding to features) are the rows (but should be columns)

    scaling: - 'standard': use standard scaler to scale the features
             - 'log2': use log2 eof all values
             - 'log2-standard' First apply log2, then standardize
             - None: don't scale
             'standard' and 'log2' help to avoid PCA domination by large genes

    '''
    if dimension not in [2]:  # for now only support 2 dimension plot
        raise ValueError('Only 2D plot supported')

    pca_data = df.T

    if 'log2' in scaling:
        pca_data = np.log2(pca_data + 1)

    if 'standard' in scaling:
        scaler = StandardScaler()
        pca_data = scaler.fit_transform(pca_data)

    # df.index.name = 'miRNA'
    pca = PCA()
    pca.fit(pca_data)

    explained_variance = np.sum(pca.explained_variance_ratio_[:dimension])

    ax = plot_pca_2d_projection(pca, pca_data, df.columns, ax=ax, cmap=cmap, label_dots=label_dots)

    ax.set_title(f'Explained variance: {explained_variance}')
    return ax

df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0)

fig, ax = plt.subplots(figsize=(16, 16))
pcaplot(df, scaling='log2', ax=ax)

fig.savefig(snakemake.output[0])
