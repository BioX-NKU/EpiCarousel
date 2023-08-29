import sys
import numpy as np
#from memory_profiler import profile
import igraph as ig
import anndata as ad
# from scipy import sparse
from scipy.sparse import csr_matrix, vstack
import preprocess as pp
import warnings
import gc
warnings.filterwarnings("ignore")
gc.collect()

def identify_metacells(data_name, 
                          K: int,
                          chunk_dir: str,    
                          chunk_preprocessed_dir:str,    
                          chunk_mc_dir: str, 
                          chunk_mc_binarized_dir: str,
                          carousel_resolution: int = 10,
                          steps: int = 4,
                          if_bi: int = 1,
                          if_mc_bi: int = 1,
                          mc_mode: str = 'average',
                          threshold: float = 0.0,
                          filter_rate: float = 0.01,
                          neighbors_method: str = 'umap',
                          n_components: int = 50, 
                          svd_solver: str = 'arpack'
                         ):
    """
    Identify metacells.
    """

    adata = pp.read(chunk_dir + "/%s_fold%d.h5ad"%(data_name,K))

    adataX_origin = adata.X.copy()
    adata = pp.ATAC_preprocess(adata,
                               filter_rate=filter_rate,
                               if_bi=if_bi, 
                               method=neighbors_method,
                               n_components=n_components,
                               svd_solver=svd_solver
                               )
    print(adata)
    adata.write(chunk_preprocessed_dir + '/%s_fold%d_if_bi_%d.h5ad'%(data_name, K, if_bi))

    # Graph constrution
    graph = ig.Graph.Weighted_Adjacency(adata.obsp['connectivities'],mode = 'undirected')


    # Walktrap
    community = graph.community_walktrap(steps = steps, weights = graph.es.get_attribute_values('weight'))
    communities = community.as_clustering(n = int(adata.X.shape[0]/carousel_resolution))

    del adata
    gc.collect()
    
    # Metacell identification 
    out = []
    if mc_mode == 'average':
        for group, idx in enumerate(communities):
            lim = csr_matrix(np.ravel((np.sum(adataX_origin[idx], axis=0) / len(idx))))
            out.append(lim)
    elif mc_mode == 'sum':
        for group, idx in enumerate(communities):
            lim = csr_matrix(np.ravel((np.sum(adataX_origin[idx], axis=0))))
            out.append(lim)
        
    mc_adata = ad.AnnData(X=vstack(out).astype(np.float32))
    mc_adata.X = csr_matrix(mc_adata.X)
    
    dwt = []
    for i in range(mc_adata.X.shape[0]):
        dwt.append(K)
    mc_adata.obs['which_fold'] = dwt
    dwt = []
    #mc_adata_cp.obs['cells'] = dwt
    for i, community in enumerate(communities):
        dwt.append(str(community).replace(',','').replace('[','').replace(']',''))
    mc_adata.obs['cells'] = dwt
    if if_mc_bi == 0:    
        mc_adata.write(chunk_mc_dir + "/%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad"%(data_name,K, if_bi, if_mc_bi))

    if if_mc_bi == 1:
        if threshold == 0:
            mc_adata.X.data = np.ones(mc_adata.X.data.shape[0], dtype=np.int8)
        else:
            from sklearn.preprocessing import Binarizer
            binarizer = Binarizer(threshold=threshold)
            mc_adata.X = binarizer.transform(mc_adata.X).astype(np.int8)
        mc_adata.write(chunk_mc_binarized_dir + "/%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (data_name, K, if_bi, if_mc_bi))

    del mc_adata
    gc.collect()

    
if __name__ == '__main__':
    global K
    data_name = sys.argv[1]
    K = np.int(sys.argv[2])
    chunk_dir = sys.argv[3]
    chunk_preprocessed_dir = sys.argv[4]
    chunk_mc_dir = sys.argv[5]
    chunk_mc_binarized_dir = sys.argv[6]
    carousel_resolution = np.int(sys.argv[7])
    steps = np.int(sys.argv[8])
    if_bi = np.int(sys.argv[9])
    if_mc_bi = np.int(sys.argv[10])
    mc_mode = sys.argv[11]
    threshold = np.float(sys.argv[12])
    filter_rate = np.float(sys.argv[13])
    neighbors_method = sys.argv[14]
    n_components = np.int(sys.argv[15])
    svd_solver = sys.argv[16]
    
    identify_metacells(data_name, K, chunk_dir, chunk_preprocessed_dir, chunk_mc_dir, chunk_mc_binarized_dir, carousel_resolution, steps, if_bi, if_mc_bi,mc_mode, threshold, filter_rate, neighbors_method, n_components, svd_solver)
    