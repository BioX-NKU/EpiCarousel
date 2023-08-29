import scanpy as sc
import episcanpy.api as epi
import numpy as np
from scipy.sparse import csr_matrix


def read(filename,
         formal: str = 'normal'):
    """
    Load the h5ad file.
    """    
    if formal == 'normal':
        adata = sc.read(filename)
    elif formal == 'lazily':
        import snapatac2 as snap
        adata = snap.read(filename, backed='r')
    else:
        raise RuntimeError('Formal %s is invalid!' % formal)
    return adata

from sklearn.feature_extraction.text import TfidfTransformer
def tfidf(count_mat): 
    """
    Perform TF-IDF transformation.
    """    
    model = TfidfTransformer(smooth_idf=False, norm="l2")
    model = model.fit(np.transpose(count_mat))
    model.idf_ -= 1
    tf_idf = np.transpose(model.transform(np.transpose(count_mat))).astype(np.float32)
    return csr_matrix(tf_idf)

def ATAC_preprocess(adata,
                    filter_rate: float = 0.01,
                    n: int = 15,
                    metric: str = 'euclidean',
                    method='umap',
                    if_bi: int = 0,
                    n_components: int = 50, 
                    svd_solver: str = 'arpack'
      ):
    """
    Preprocess scCAS data matrix.

    Parameters
    ----------
    adata :  AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to peaks/regions.
    filter_rate : float, optional
        Proportion for feature selection, by default 0.01
    """
    if if_bi == 1:
        adata.X.data = np.ones(adata.X.data.shape[0], dtype = np.int8)

    epi.pp.filter_features(adata, min_cells=np.ceil(filter_rate*adata.shape[0]))
    
    # TF-IDF transformation
    adata.X = csr_matrix(tfidf(adata.X.T).T )
    epi.pp.pca(adata, svd_solver=svd_solver, n_comps=n_components)
    epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method, use_rep='X_pca')

    return adata
