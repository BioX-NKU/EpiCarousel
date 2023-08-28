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

def metacell_construction(data_name, 
                          K: int,
                          chunk_dir: str,    # chunk_dir
                          chunk_preprocessed_dir:str,    # chunk_preprocessed_dir
                          chunk_mc_dir: str,  # mc_adata_K保存地址
                          chunk_mc_binarized_dir: str, # 二值化后的mc_adata_K保存地址
                          chunk_method: str = 'tfidf3',
                          carousel_resolution: int = 10,
                          steps: int = 4,
                          if_bi: int = 1,
                          if_mc_bi: int = 1,
                          mc_mode: str = 'average',
                          threshold: float = 0.0,
                          filter_rate: float = 0.01,
                          neighbors_method: str = 'umap',
                          decomposition: str = 'pca',
                          n_components: int = 50, 
                          svd_solver: str = 'arpack',
                          random_state: int = 1,
                          n_iter: int = 7,
                          kernel: str = 'linear'
                         ):

    adata = pp.read(chunk_dir + "/%s_fold%d.h5ad"%(data_name,K))

    adataX_origin = adata.X.copy()
    adata = pp.ATAC_preprocess(adata, 
                               filter_rate=filter_rate, 
                               transform=chunk_method, 
                               if_bi=if_bi, 
                               method=neighbors_method, 
                               decomposition=decomposition,
                              n_components=n_components,
                              svd_solver=svd_solver,
                              random_state=random_state,
                              n_iter=n_iter,
                              kernel=kernel)
    print(adata)
    adata.write(chunk_preprocessed_dir + '/%s_fold%d_%s_if_bi_%d.h5ad'%(data_name, K, chunk_method, if_bi))

    # graph constrution
    graph = ig.Graph.Weighted_Adjacency(adata.obsp['connectivities'],mode = 'undirected')


    # walktrap
    community = graph.community_walktrap(steps = steps, weights = graph.es.get_attribute_values('weight'))
    communities = community.as_clustering(n = int(adata.X.shape[0]/carousel_resolution))

    del adata
    gc.collect()
    
    #metacell identification 
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
        mc_adata.write(chunk_mc_dir + "/%s_chunk_%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad"%(data_name, chunk_method,K, if_bi, if_mc_bi))

    if if_mc_bi == 1:
        if threshold == 0:
            mc_adata.X.data = np.ones(mc_adata.X.data.shape[0], dtype=np.int8)
        else:
            from sklearn.preprocessing import Binarizer
            binarizer = Binarizer(threshold=threshold)
            mc_adata.X = binarizer.transform(mc_adata.X).astype(np.int8)
        mc_adata.write(chunk_mc_binarized_dir + "/%s_chunk_%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (data_name, chunk_method, K, if_bi, if_mc_bi))

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
    chunk_method = sys.argv[7]
    carousel_resolution = np.int(sys.argv[8])
    steps = np.int(sys.argv[9])
    if_bi = np.int(sys.argv[10])
    if_mc_bi = np.int(sys.argv[11])
    mc_mode = sys.argv[12]
    threshold = np.float(sys.argv[13])
    filter_rate = np.float(sys.argv[14])
    neighbors_method = sys.argv[15]
    decomposition = sys.argv[16]
    n_components = np.int(sys.argv[17])
    svd_solver = sys.argv[18]
    random_state = np.int(sys.argv[19])
    n_iter = np.int(sys.argv[20])
    kernel = sys.argv[21]
    
    metacell_construction(data_name, K, chunk_dir, chunk_preprocessed_dir, chunk_mc_dir, chunk_mc_binarized_dir, chunk_method, carousel_resolution, steps, if_bi, if_mc_bi,mc_mode, threshold, filter_rate, neighbors_method, decomposition, n_components, svd_solver, random_state, n_iter, kernel)
    