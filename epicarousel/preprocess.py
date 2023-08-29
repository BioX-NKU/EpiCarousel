import scanpy as sc
import episcanpy.api as epi
import numpy as np
import gc
import snapatac2 as snap
from scipy.sparse import csr_matrix


def read(filename,
         formal: str = 'normal'):
    
    if formal == 'normal':
        adata = sc.read(filename)
    elif formal == 'lazily':
        adata = snap.read(filename, backed='r')
    else:
        raise RuntimeError('Formal %s is invalid!' % formal)
    return adata

from sklearn.feature_extraction.text import TfidfTransformer
def tfidf3(count_mat): 
    model = TfidfTransformer(smooth_idf=False, norm="l2")
    model = model.fit(np.transpose(count_mat))
    model.idf_ -= 1
    tf_idf = np.transpose(model.transform(np.transpose(count_mat))).astype(np.float32)
    return csr_matrix(tf_idf)

def ATAC_preprocess(adata,
                    filter_rate=0.01,
                    transform: str = 'tfidf3',
                    n: int = 15,
                    metric: str = 'euclidean',
                    method='umap',
                    if_bi: int = 0,
                    decomposition: str = 'pca',
                    n_components: int = 50, 
                    svd_solver: str = 'arpack',
                    random_state: int = 1,
                    n_iter: int = 7,
                    kernel: str = 'linear'
      ):


    if if_bi == 1:
        adata.X.data = np.ones(adata.X.data.shape[0], dtype = np.int8)

    epi.pp.filter_features(adata, min_cells=np.ceil(filter_rate*adata.shape[0]))
    
    if transform == 'tfidf1':
        adata.X = tfidf1(adata.X.T).T 
#         epi.pp.pca(adata, svd_solver='arpack')
#         epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method)
    elif transform == 'tfidf2':
        adata.X = tfidf2(adata.X.T).T 
#         epi.pp.pca(adata, svd_solver='arpack')
#         epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method)
    elif transform == 'tfidf3' or transform == 'tfidf':
        adata.X = csr_matrix(tfidf3(adata.X.T).T )
#         epi.pp.pca(adata, svd_solver='arpack')
#         epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method)
    elif transform == 'tfidf22':
        adata = tfidf22(adata)
#         epi.pp.pca(adata, svd_solver='arpack')
#         epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method)
    elif transform == 'normalize_total_log_transform':
        sc.pp.normalize_total(adata)
        epi.pp.log1p(adata)
#         epi.pp.pca(adata, svd_solver='arpack')
#         epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method)
    elif transform == 'normalize_total':
        sc.pp.normalize_total(adata)
#         epi.pp.pca(adata, svd_solver='arpack')
#         epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method)
    else:
        raise RuntimeError('Cannot transform by '+ transform + '!')
        
    if decomposition == 'pca':
        epi.pp.pca(adata, svd_solver=svd_solver, n_comps=n_components)
        epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method, use_rep='X_pca')
    elif decomposition == 'LSA':
        print(" Start LSA!")
        from sklearn.decomposition import TruncatedSVD
        svd = TruncatedSVD(n_components=n_components, n_iter=n_iter, random_state=random_state)
        adata.obsm['X_svd'] = svd.fit_transform(adata.X)
        del svd
        gc.collect()
        epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method, use_rep='X_svd')
    elif decomposition == 'KernelPCA':
        from sklearn.decomposition import KernelPCA
        transformer = KernelPCA(n_components=n_components, kernel=kernel, random_state=random_state)
        adata.obsm['X_KernelPCA'] = transformer.fit_transform(adata.X)
        del transformer
        gc.collect()
        epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method, use_rep='X_KernelPCA')
    elif decomposition == 'IterativeLSI':
        print("Start IterativeLSI.")
        from gensim.models import LsiModel
        def LSI(ad, n_topics = 50):
            corpus = []
            for i in range(len(ad.X.indptr)-1):
                corpus_i = []
                for j in range(ad.X.indptr[i],ad.X.indptr[i+1]):
                    corpus_i.append((ad.X.indices[j],ad.X.data[j]))
                corpus.append(corpus_i)


            model = LsiModel(corpus,num_topics=n_topics)
            vectorized_corpus = model[corpus]
            nparray = np.zeros((ad.shape[0],n_topics))
            assert ad.shape[0] == len(vectorized_corpus)
            k=-1
            for vector in vectorized_corpus:
                k = k+1
                for topic in vector:
                    nparray[k][topic[0]] = topic[1]
            ad.obsm['X_LSI'] = nparray
            del model
            del vectorized_corpus
            del nparray
            del corpus
            gc.collect()
            return ad
        LSI(adata, n_topics=n_components)
        epi.pp.neighbors(adata, n_neighbors=n, metric=metric, method=method, use_rep='X_LSI')
    return adata
