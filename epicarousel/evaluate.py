import pandas as pd
import numpy as np

def get_purity(obs, mc_adata, index, chunk_size):
    """
    Calculate the purity of each metacell.
    """    
    purity_list = []
    for i in range(mc_adata.shape[0]):
        cell_list = mc_adata.obs['cells'][i].split()
        fold = int(mc_adata.obs['which_fold'][i])
        idx = [(int(item) + (fold-1) * chunk_size) for item in cell_list]
        purity_list.append(pd.value_counts(obs[index][idx]).max()/pd.value_counts(obs[index][idx]).sum())
    return purity_list

def get_celltype(obs, mc_adata, index, chunk_size):
    """
    Get the cell types of single cells contained in corresponding metacells.
    """    
    celltype_list = []
    for i in range(mc_adata.shape[0]):
        cell_list = mc_adata.obs['cells'][i].split()
        fold = int(mc_adata.obs['which_fold'][i])
        idx = [(int(item) + (fold-1) * chunk_size) for item in cell_list]
        celltype_list.append(pd.value_counts(obs[index][idx]).index.tolist()[0])
    return celltype_list

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
    
    print('AMI:%.3f\tARI:%.3f\tNMI:%.3f\tHomo:%.3f\tCS:%.3f\tV-measure:%.3f\tFMI:%.3f\t'%(AMI,ARI,NMI,HOM,completeness,VMS,FMS))
    return AMI, ARI, NMI, HOM, completeness, VMS, FMS


