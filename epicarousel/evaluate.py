import pandas as pd
import numpy as np

def get_purity(obs, mc_adata, index, chunk_size):
    purity_list = []
    for i in range(mc_adata.shape[0]):
        cell_list = mc_adata.obs['cells'][i].split()
        fold = int(mc_adata.obs['which_fold'][i])
        idx = [(int(item) + (fold-1) * chunk_size) for item in cell_list]
        purity_list.append(pd.value_counts(obs[index][idx]).max()/pd.value_counts(obs[index][idx]).sum())
    return purity_list

def get_celltype(obs, mc_adata, index, chunk_size):
    purity_list = []
    for i in range(mc_adata.shape[0]):
        cell_list = mc_adata.obs['cells'][i].split()
        fold = int(mc_adata.obs['which_fold'][i])
        idx = [(int(item) + (fold-1) * chunk_size) for item in cell_list]
        purity_list.append(pd.value_counts(obs[index][idx]).index.tolist()[0])
    return purity_list

def cluster_evaluation(adata_obs, label_key, cluster_key):
    '''
    Clustering Performance Evaluation

    Args:
        adata_obs: polars.internals.frame.DataFrame.
        label_key: e.g. 'cell type', 'cell_type'
        cluster_key: e.g. 'mc_Dleiden'

    Returns:
        evaluation

    '''
    from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, normalized_mutual_info_score, homogeneity_score, completeness_score, v_measure_score, fowlkes_mallows_score
    ARI = adjusted_rand_score(adata_obs[label_key], adata_obs[cluster_key])  # ARI is a symmetric measure
    AMI = adjusted_mutual_info_score(adata_obs[label_key], adata_obs[cluster_key])  # This metric is furthermore symmetric
    NMI = normalized_mutual_info_score(adata_obs[label_key], adata_obs[cluster_key])  # This metric is furthermore symmetric
    HOM = homogeneity_score(adata_obs[label_key], adata_obs[cluster_key])  # This metric is not symmetric
    completeness = completeness_score(adata_obs[label_key], adata_obs[cluster_key])  # This metric is not symmetric
    VMS = v_measure_score(adata_obs[label_key], adata_obs[cluster_key])  # This metric is furthermore symmetric
    FMS = fowlkes_mallows_score(adata_obs[label_key], adata_obs[cluster_key])
    
    print('AMI:%.3f\tARI:%.3f\tNMI:%.3f\tHOM:%.3f\tcompleteness:%.3f\tV-measure:%.3f\tFMS:%.3f\t'%(AMI,ARI,NMI,HOM,completeness,VMS,FMS))
    return AMI, ARI, NMI, HOM, completeness, VMS, FMS

def silhouette_score(X,
                     labels,
                     metric='euclidean',
                     sample_size: int = None,
                     random_state: int = 1):

    from sklearn.metrics import silhouette_score
    return(silhouette_score(X, labels=labels, metric=metric, sample_size=sample_size, random_state=random_state))

def compactness(ad, low_dim_embedding='X_pca', metacell_label='carousel'):
    """
    Compute compactness of each metacell. Compactness is defined as the average variance of diffusion components
    across cells that constitute a metcell

    :param ad: (Anndata) Anndata object, source_adata, not metacell_adata
    :param low_dim_embedding: (str) `ad.obsm` field for constructing diffusion components
    :param metacell_label: (str) `ad.obs` field for computing diffusion component variances

    :return: `pd.DataFrame` with a dataframe of compactness per metacell

    """

    import palantir

    components = pd.DataFrame(ad.obsm[low_dim_embedding]).set_index(ad.obs_names)
    dm_res = palantir.utils.run_diffusion_maps(components)
    dc = palantir.utils.determine_multiscale_space(dm_res, n_eigs=10)

    return pd.DataFrame(dc.join(ad.obs[metacell_label]).groupby(metacell_label).var().mean(1)).rename(columns={0:'compactness'})

def separation(ad,
                                   low_dim_embedding='X_pca',
                                   nth_nbr=1,
                                   cluster=None,
                                   carousels_label='carousel'):
    """
    Compute separation of each metacell. Separation is defined is the distance to the nearest neighboring metacell

    :param ad: (Anndata) Anndata object, source_adata, not metacell_adata
    :param low_dim_embedding: (str) `ad.obsm` field for constructing diffusion components
    :param nth_nbr: (int) Which neighbor to use for computing separation
    :param carousel_label: (str) `ad.obs` field for computing diffusion component variances

    :return: `pd.DataFrame` with a separation of compactness per metacell

    """
    import palantir

    components = pd.DataFrame(ad.obsm[low_dim_embedding]).set_index(ad.obs_names)
    dm_res = palantir.utils.run_diffusion_maps(components)
    dc = palantir.utils.determine_multiscale_space(dm_res, n_eigs=10)

    # Compute DC per metacell
    metacells_dcs = dc.join(ad.obs[carousels_label], how='inner').groupby(carousels_label).mean()

    from sklearn.neighbors import NearestNeighbors
    neigh = NearestNeighbors(n_neighbors=nth_nbr)
    nbrs = neigh.fit(metacells_dcs)
    dists, nbrs = nbrs.kneighbors()
    dists = pd.DataFrame(dists).set_index(metacells_dcs.index)
    dists.columns += 1

    nbr_cells = np.array(metacells_dcs.index)[nbrs]

    metacells_nbrs = pd.DataFrame(nbr_cells)
    metacells_nbrs.index = metacells_dcs.index
    metacells_nbrs.columns += 1

    if cluster is not None:

        # Get cluster type of neighbors to ensure they match the metacell cluster
        clusters = ad.obs.groupby(carousels_label)[cluster].agg(lambda x: x.value_counts().index[0])
        nbr_clusters = pd.DataFrame(clusters.values[nbrs]).set_index(clusters.index)
        nbr_clusters.columns = metacells_nbrs.columns
        nbr_clusters = nbr_clusters.join(pd.DataFrame(clusters))

        clusters_match = nbr_clusters.eq(nbr_clusters[cluster], axis=0)
        return pd.DataFrame(dists[nth_nbr][clusters_match[nth_nbr]]).rename(columns={1:'separation'})
    else:
        return pd.DataFrame(dists[nth_nbr]).rename(columns={1:'separation'})

