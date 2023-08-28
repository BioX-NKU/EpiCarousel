import episcanpy.api as epi
import pandas as pd
import numpy as np
import random
import anndata as ad
from tqdm import tqdm
import os
from shutil import rmtree
from sklearn.utils import shuffle
from memory_profiler import profile
from subprocess import run
from pkg_resources import resource_filename

from . import preprocess as pp
from . import evaluate

import warnings
import gc

warnings.filterwarnings("ignore")
gc.collect()


def setup_seed(seed):
    """
    Set random seed.

    Parameters
    ----------
    seed
        Number to be set as random seed for reproducibility.

    """
    np.random.seed(seed)
    random.seed(seed)


# setup_seed(1)

def constrcution(fold_number: int,
                 data_name: str,
                 chunk_dir: str,
                 chunk_preprocessed_dir: str,
                 chunk_mc_dir: str,
                 chunk_mc_binarized_dir: str,
                 chunk_method: str,
                 carousel_resolution: int,
                 step: int,
                 threads: int,
                 if_bi: int,
                 if_mc_bi: int,
                 mc_mode: str,
                 threshold: float,
                 filter_rate: float,
                 neighbors_method: str,
                 decomposition: str,
                 n_components: int,
                 svd_solver: str,
                 random_state: int,
                 n_iter: int,
                 kernel: str
                 ):
    """
    Identify metacells in scCAS data by EpiCarousel.

    Parameters
    ----------
    fold_number : int
        The ordering of the input chunk data in initial single cells data.
    data_name : str
        The name of the datasets.
    chunk_dir : str
        Location of the saved chunk files.
    chunk_preprocessed_dir : str
        Location of the preprocessed chunk files.
    chunk_mc_dir : str
        Location of metacells h5ad file storage for the chunk.
    chunk_mc_binarized_dir : str
        Location of the binarized metacells h5ad file storage for the chunk.
    chunk_method : str
        The version of TF-IDF transformation to be performed on chunk matrix.
    carousel_resolution : int
        The ratio of the number of original cells to that of metacells.
    step : int
        Walktrap length.
    threads : int
        _description_
    if_bi : int
        _description_
    if_mc_bi : int
        _description_
    mc_mode : str
        _description_
    threshold : float
        _description_
    filter_rate : float
        _description_
    neighbors_method : str
        _description_
    decomposition : str
        _description_
    n_components : int
        _description_
    svd_solver : str
        _description_
    random_state : int
        _description_
    n_iter : int
        _description_
    kernel : str
        _description_
    """    """"""        
    package_path = resource_filename('epicarousel', "")
    os.chmod(package_path + '/construction.sh', 0o755)
    print(os.getcwd())
    os.chdir(package_path)
    print(os.getcwd())
    print(package_path)
    os.system('./construction.sh '
              + str(fold_number)
              + ' ' + data_name
              + ' ' + chunk_dir
              + ' ' + chunk_preprocessed_dir
              + ' ' + chunk_mc_dir
              + ' ' + chunk_mc_binarized_dir
              + ' ' + chunk_method
              + ' ' + str(carousel_resolution)
              + ' ' + str(step)
              + ' ' + str(threads)
              + ' ' + str(if_bi)
              + ' ' + str(if_mc_bi)
              + ' ' + mc_mode
              + ' ' + str(threshold)
              + ' ' + str(filter_rate)
              + ' ' + neighbors_method
              + ' ' + decomposition
              + ' ' + str(n_components)
              + ' ' + svd_solver
              + ' ' + str(random_state)
              + ' ' + str(n_iter)
              + ' ' + kernel
             )


# @profile
def consolidate_labels(atac_ad,
                       atac_mc_ad,
                       chunk_size):
    # 赋mc标签
    mc_label = ['carousel']
    for s in mc_label:
        label = atac_mc_ad.obs.index.tolist()
        label.append('None')
        dwt = ['None'] * atac_ad.shape[0]
        dwt = pd.Categorical(dwt, categories=label)
        atac_ad.obs['%s' % s] = dwt
        for i in range(atac_mc_ad.shape[0]):
            fold = atac_mc_ad.obs['which_fold'][i]
            cell_list = atac_mc_ad.obs['cells'][i].split(' ')
            cell_list = [(int(item) + (fold - 1) * chunk_size) for item in cell_list]
            atac_ad.obs['%s' % s][cell_list] = atac_mc_ad.obs.index[i]
        atac_ad.obs['%s' % s] = pd.Categorical(atac_ad.obs['%s' % s],
                                               categories=list(atac_ad.obs['%s' % s].unique()))
    # atac_ad.obs['carousel']

    # 赋mc聚类标签
    cluster_name = ['Dleiden', 'Dlouvain', 'Cleiden', 'Clouvain']
    for s in cluster_name:
        cluster_label = atac_mc_ad.obs[s].cat.categories.tolist()
        cluster_label.append('None')
        dwt = ['None'] * atac_ad.shape[0]
        dwt = pd.Categorical(dwt, categories=cluster_label)
        atac_ad.obs['mc_%s' % s] = dwt
        for i in range(atac_mc_ad.shape[0]):
            fold = atac_mc_ad.obs['which_fold'][i]
            cell_list = atac_mc_ad.obs['cells'][i].split(' ')
            cell_list = [(int(item) + (fold - 1) * chunk_size) for item in cell_list]
            atac_ad.obs['mc_%s' % s][cell_list] = atac_mc_ad.obs[s][i]
        atac_ad.obs['mc_%s' % s] = pd.Categorical(atac_ad.obs['mc_%s' % s],
                                                  categories=list(atac_ad.obs['mc_%s' % s].unique()))
    # atac_ad.obs[['mc_Dleiden', 'mc_Dlouvain', 'mc_Cleiden', 'mc_Clouvain']]

    # 赋fold标签
    mc_label = ['which_fold']
    for s in mc_label:
        label = list(atac_mc_ad.obs['which_fold'].unique())
        label.append('None')
        dwt = ['None'] * atac_ad.shape[0]
        dwt = pd.Categorical(dwt, categories=label)
        atac_ad.obs[s] = dwt
        for i in range(atac_mc_ad.shape[0]):
            fold = atac_mc_ad.obs['which_fold'][i]
            cell_list = atac_mc_ad.obs['cells'][i].split(' ')
            cell_list = [(int(item) + (fold - 1) * chunk_size) for item in cell_list]
            atac_ad.obs[s][cell_list] = atac_mc_ad.obs['which_fold'][i]
        atac_ad.obs[s] = pd.Categorical(atac_ad.obs[s], categories=list(atac_ad.obs[s].unique()))
    # atac_ad.obs['which_fold']

    return atac_ad


class Carousel:

    def __init__(self,
                 data_name: str,
                 data_dir: str,
                 if_bi: int = 1,
                 if_mc_bi: int = 1,
                 threshold: float = 0.0,
                 filter_cells: bool = False,
                 filter_rate: float = 0.01,
                 chunk_size: int = 10000,
                 original_method: str = 'tfidf3',
                 chunk_method: str = 'tfidf3',
                 carousel_method: str = 'tfidf3',
                 carousel_resolution: int = 10,
                 base: str = '/home/metacell/data/metacell/carousel/output',
                 shuffle: int = 0,
                 random_state: int = 1,
                 step: int = 4,
                 threads: int = 8,
                 index: str = 'cell_type',
                 mc_mode: str = 'average',
                 neighbors_method: str = 'umap',
                 decomposition: str = 'pca',
                 n_components: int = 50,
                 svd_solver: str = 'arpack',
                 n_iter: int = 7,
                 kernel: str = 'linear'
                 ):
        '''

        Args:
            data_name: The name of data.
            data_dir: The path where the data is saved.
            if_bi: int, whether the original data is binarized, 1 represents binarization
                and other numbers represent no binarization
            if_mc_bi: int, whether metacell count matrix is binared, 1 represents binarization
                        and other numbers represent no binarization
            threshold: float, the threshold for metacell count matrix
            filter_cells:
            filter_rate: min_cells=np.ceil(filter_rate*adata.shape[0])
            chunk_size:
            original_method:
            chunk_method:
            carousel_method:
            carousel_resolution:
            base:
            shuffle:
            random_state:
            step:
            threads:
            index:
            mc_mode:
            neighbors_method:
        '''
        self.data_name = data_name
        self.filter_cells = filter_cells
        self.chunk_size = chunk_size
        self.original_method = original_method  # 原数据直接聚类的预处理方法
        self.chunk_method = chunk_method
        self.metacell_method = carousel_method
        self.carousel_resolution = carousel_resolution
        self.gamma = 1 / carousel_resolution
        self.if_bi = if_bi  # 原数据是否做二值化
        self.if_mc_bi = if_mc_bi  # 控制输出的metacell数据是否二值化过
        self.threshold = threshold
        self.random_state = random_state
        self.step = step
        self.threads = threads  # 并行数
        self.shuffle = shuffle
        self.index = index
        self.mc_mode = mc_mode
        self.filter_rate = filter_rate
        self.neighbors_method = neighbors_method
        self.decomposition = decomposition
        self.n_components = n_components
        self.svd_solver = svd_solver
        self.n_iter = n_iter
        self.kernel = kernel

        self.adata = None
        self.cells_number = None
        self.peaks_number = None
        self.fold_number = None
        self.mc_adata = None
        self.sparsity = None
        self.mc_sparsity = None

        # dir
        self.base = base
        self.data_dir = data_dir
        self.proj_dir = None
        self.sf_dir = None
        self.chunk_dir = None
        self.chunk_preprocessed_dir = None
        self.chunk_mc_dir = None
        self.chunk_mc_binarized_dir = None
        self.carousel_dir = None
        self.carousel_binarized_dir = None
        self.carousel_preprocessed_dir = None
        # self.carousel_filtered_dir = None
        self.carousel_clustered_dir = None
        self.result_dir = None
        self.dirs = None

        self.var = None
        self.obs = None
        self.result = None
        self.res = None

#     @profile
    def make_dirs(self):
        '''


        '''
        self.proj_dir = self.base + '/' + self.data_name + '/' + str(self.chunk_size) + '_is_bi_' + str(
            self.if_bi) + '_' + self.chunk_method + '_' + '_steps' + str(self.step) + '_resolution_' + str(
            self.carousel_resolution) + '_mc_mode_' + self.mc_mode + '_if_mc_bi_' + str(
            self.if_mc_bi) + '_threshold_' + str(
            self.threshold) + '_' + self.metacell_method + '_' + self.neighbors_method + '_export'
        print("Project output directory: %s" % self.proj_dir)
        self.sf_dir = self.proj_dir + '/shuffled'
        self.chunk_dir = self.proj_dir + '/chunk'
        self.chunk_preprocessed_dir = self.proj_dir + '/chunk_preprocessed'
        self.chunk_mc_dir = self.proj_dir + '/chunk_mc'
        self.chunk_mc_binarized_dir = self.proj_dir + '/chunk_mc_binarized'
        self.carousel_dir = self.proj_dir + '/carousel'
        self.carousel_binarized_dir = self.proj_dir + '/carousel_binarized'
        self.carousel_preprocessed_dir = self.proj_dir + '/carousel_pped'
        self.carousel_clustered_dir = self.proj_dir + '/carousel_clustered'
        self.result_dir = self.proj_dir + '/carousel_clustering_performance'

        self.dirs = [self.proj_dir, self.sf_dir, self.chunk_dir, self.chunk_preprocessed_dir, self.chunk_mc_dir,
                     self.chunk_mc_binarized_dir, self.carousel_dir, self.carousel_binarized_dir,
                     self.carousel_preprocessed_dir, self.carousel_clustered_dir, self.result_dir]
        for directory in self.dirs:
            if not os.path.isdir(directory):
                os.makedirs(directory)

#     @profile
    def shuffle_data(self):

        print("Start shuffling data!")
        if not os.path.exists(self.sf_dir + '/%s_sf_%d.h5ad' % (self.data_name, self.random_state)):
            adata = pp.read(self.data_dir)
            adata = shuffle(adata, random_state=self.random_state)
            adata.write(self.sf_dir + '/%s_sf_%d.h5ad' % (self.data_name, self.random_state))
            del adata
        gc.collect()
        print("Finish shuffling data!")

#     @profile
    def data_split(self):

        print("Start data_split!")
        if self.shuffle == 1:
            snap_adata = pp.read(self.sf_dir + '/%s_sf_%s.h5ad' % (self.data_name, str(self.random_state)),
                                 formal='lazily')
        else:
            snap_adata = pp.read(self.data_dir, formal='lazily')

        snap_adata.X.enable_cache()

        chunks = snap_adata.X.chunked(self.chunk_size)
        self.cells_number, self.peaks_number = snap_adata.shape
        self.obs = snap_adata.obs[:]
        self.obs = pd.DataFrame(self.obs[:, 0:], index=self.obs.columns[0:]).T
        self.var = snap_adata.var[:]
        self.var = pd.DataFrame(self.var[:, 0:], index=self.var.columns[0:]).T
        origin_sum = 0

        for i, chunk in tqdm(enumerate(chunks)):
#             chunk_adata = ad.AnnData(chunk[0])
#             chunk_adata.X = chunk_adata.X.astype(int)
            chunk_adata = ad.AnnData(X=chunk[0].astype(np.int8))
            chunk_adata.write(self.chunk_dir + "/%s_fold%d.h5ad" % (self.data_name, i + 1))
            origin_sum += np.sum(chunk_adata.X)
            del chunk_adata
            del chunk
            gc.collect()

        self.fold_number = int(np.ceil(self.cells_number / self.chunk_size))
        self.sparsity = origin_sum / self.cells_number / self.peaks_number
        # self.chunk_size = chunk_size
        snap_adata.close()
        del origin_sum
        gc.collect()
        print("Finish data_split!")

#     @profile
    def original_data_read(self):

        if self.shuffle:
            self.adata = pp.read(self.sf_dir + '/%s_sf_%d.h5ad' % (self.data_name, self.random_state))
        else:
            self.adata = pp.read(self.data_dir)

    @profile
    def original_data_preprocess(self, method='umap'):

        self.adata = pp.ATAC_preprocess(self.adata,
                                        filter_rate=0.01,
                                        transform=self.original_method,
                                        n=15,
                                        metric='euclidean',
                                        method=method,
                                        if_bi=0
                                        )

#     @profile
    def original_data_clustering(self):

        epi.tl.leiden(self.adata, key_added='Dleiden')

        epi.tl.louvain(self.adata, key_added='Dlouvain')

        epi.tl.getNClusters(self.adata, n_cluster=self.adata.obs[self.index].nunique(), method='leiden')
        self.adata.obs['Cleiden'] = self.adata.obs['leiden'].copy()

        epi.tl.getNClusters(self.adata, n_cluster=self.adata.obs[self.index].nunique(), method='louvain')
        self.adata.obs['Clouvain'] = self.adata.obs['louvain'].copy()

        # del self.adata
        gc.collect()

#     @profile
    def metacell_construction(self):

        constrcution(fold_number=self.fold_number,
                     data_name=self.data_name,
                     chunk_dir=self.chunk_dir,
                     chunk_preprocessed_dir=self.chunk_preprocessed_dir,
                     chunk_mc_dir=self.chunk_mc_dir,
                     chunk_mc_binarized_dir=self.chunk_mc_binarized_dir,
                     chunk_method=self.chunk_method,
                     carousel_resolution=self.carousel_resolution,
                     step=self.step,
                     threads=self.threads,
                     if_bi=self.if_bi,
                     if_mc_bi=self.if_mc_bi,
                     mc_mode=self.mc_mode,
                     threshold=self.threshold,
                     filter_rate=self.filter_rate,
                     neighbors_method=self.neighbors_method,
                     decomposition=self.decomposition,
                     n_components=self.n_components,
                     svd_solver=self.svd_solver,
                     random_state=self.random_state,
                     n_iter=self.n_iter,
                     kernel=self.kernel
                     )

#     @profile
    def merge_metacell(self):
        '''
        merge metacell chunks

        '''

        print("Start merging metacells.")
        if self.if_mc_bi == 1:
            i = 0
            self.mc_adata = pp.read(self.chunk_mc_binarized_dir + "/%s_chunk_%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (
            self.data_name, self.chunk_method, i + 1, self.if_bi, self.if_mc_bi))
            for i in tqdm(range(1, self.fold_number)):
                mc_adata = pp.read(self.chunk_mc_binarized_dir + "/%s_chunk_%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (
                self.data_name, self.chunk_method, i + 1, self.if_bi, self.if_mc_bi))
                self.mc_adata = ad.concat([self.mc_adata, mc_adata])
                del mc_adata
                gc.collect()

            self.mc_adata.obs_names_make_unique()
            print("Finish merging metacells.")
            try:
                for varname in self.var[:].columns:
                    self.mc_adata.var[varname] = pd.Categorical(self.var[varname])
                del self.mc_adata.var['_index']
            except:
                pass

            # Save binarized metacell data.
            self.mc_adata.write(self.carousel_binarized_dir + '/%s_carousel_binarized.h5ad' % self.data_name)
        else:
            i = 0
            self.mc_adata = pp.read(self.chunk_mc_dir + "/%s_chunk_%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (
            self.data_name, self.chunk_method, i + 1, self.if_bi, self.if_mc_bi))
            for i in tqdm(range(1, self.fold_number)):
                mc_adata = pp.read(self.chunk_mc_dir + "/%s_chunk_%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (
                self.data_name, self.chunk_method, i + 1, self.if_bi, self.if_mc_bi))
                self.mc_adata = ad.concat([self.mc_adata, mc_adata])
                del mc_adata
                gc.collect()

            self.mc_adata.obs_names_make_unique()
            print("Finish merging metacells.")
            try:
                for varname in self.var[:].columns:
                    self.mc_adata.var[varname] = pd.Categorical(self.var[varname])
                del self.mc_adata.var['_index']
            except:
                pass

            # Save unbinarized metacell data.
            self.mc_adata.write(self.carousel_dir + '/%s_carousel.h5ad' % self.data_name)


#     @profile
    def metacell_preprocess(self,
                           neighbors_method: str = 'umap'
                           ):
        '''

        Args:
            neighbors_method:

        Returns:

        '''

        print("Start metacell preprocessing!")

        self.mc_sparsity = np.sum(self.mc_adata.X) / self.mc_adata.X.shape[0] / self.mc_adata.X.shape[1]

        self.mc_adata = pp.ATAC_preprocess(self.mc_adata,
                                           transform=self.metacell_method,
                                           filter_rate=0.01 * (self.mc_sparsity / self.sparsity),
                                           method=neighbors_method,
                                          decomposition=self.decomposition,
                                          n_components=self.n_components,
                                          svd_solver=self.svd_solver,
                                          random_state=self.random_state,
                                          n_iter=self.n_iter,
                                          kernel=self.kernel)

#         self.mc_adata.write(self.carousel_preprocessed_dir + '/%s_carousel_pped.h5ad' % self.data_name)
        print(self.mc_adata)
        print("Finish metacell preprocessing!")

#     @profile
    def metacell_data_clustering(self):

        print("Start metacell clustering!")
        print(self.obs)
        epi.tl.leiden(self.mc_adata, key_added='Dleiden')
        epi.tl.louvain(self.mc_adata, key_added='Dlouvain')

        epi.tl.getNClusters(self.mc_adata, n_cluster=len(self.obs[self.index].unique()), method='leiden')
        self.mc_adata.obs['Cleiden'] = self.mc_adata.obs['leiden'].copy()

        epi.tl.getNClusters(self.mc_adata, n_cluster=len(self.obs[self.index].unique()), method='louvain')
        self.mc_adata.obs['Clouvain'] = self.mc_adata.obs['louvain'].copy()
        print("Finish metacell clustering!")
#         self.mc_adata.write(self.carousel_clustered_dir + '/%s_carousel_clustered.h5ad' % self.data_name)

#     @profile
    def delete_dirs(self):
        '''
        删去文件和文件夹
        '''
        for direc in [self.chunk_dir, self.chunk_preprocessed_dir, self.chunk_mc_dir, self.chunk_mc_binarized_dir,
                      self.carousel_preprocessed_dir, self.carousel_clustered_dir]:
            rmtree(direc)

#     @profile
    def result_comparison(self):

        cluster_name = ['Dleiden', 'Dlouvain', 'Cleiden', 'Clouvain']
        dwt = []
        # for i in range(self.adata.X.shape[0]):
        #     dwt.append('None')

        for s in cluster_name:
            self.obs['mc_%s' % s] = ['None' for i in range(self.cells_number)]

        self.obs['carousel'] = ['None' for i in range(self.cells_number)]

        for i in range(self.mc_adata.X.shape[0]):
            cell_list = self.mc_adata.obs['cells'][i].split()
            fold = int(self.mc_adata.obs['which_fold'][i])
            cell_list = [(int(item) + (fold - 1) * self.chunk_size) for item in cell_list]
            self.obs['carousel'][cell_list] = self.mc_adata.obs_names[i]
            for s in cluster_name:
                self.obs['mc_%s' % s][cell_list] = self.mc_adata.obs['%s' % s][i]

        cluster_index = ['AMI_Dlouvain', 'ARI_Dlouvain', 'NMI_Dlouvain', 'HOM_Dlouvain', 'COM_Dlouvain', 'Vms_Dlouvain',
                         'FMS_Dlouvain',
                         'AMI_Dleiden', 'ARI_Dleiden', 'NMI_Dleiden', 'HOM_Dleiden', 'COM_Dleiden', 'Vms_Dleiden',
                         'FMS_Dleiden',
                         'AMI_Clouvain', 'ARI_Clouvain', 'NMI_Clouvain', 'HOM_Clouvain', 'COM_Clouvain', 'Vms_Clouvain',
                         'FMS_Clouvain',
                         'AMI_Cleiden', 'ARI_Cleiden', 'NMI_Cleiden', 'HOM_Cleiden', 'COM_Cleiden', 'Vms_Cleiden',
                         'FMS_Cleiden',
                         'origin_sparsity', 'carousel_sparsity', 'C-o_sparsity', 'C/o_sparsity']
        self.result = pd.DataFrame(index=range(1), columns=cluster_index, dtype=float)
        self.res = list(
            evaluate.cluster_evaluation(self.obs, self.index, 'mc_Dlouvain') + evaluate.cluster_evaluation(self.obs,
                                                                                                           self.index,
                                                                                                           'mc_Dleiden') + evaluate.cluster_evaluation(
                self.obs, self.index, 'mc_Clouvain') + evaluate.cluster_evaluation(self.obs, self.index, 'mc_Cleiden'))
        self.mc_adata.obs['purity'] = evaluate.get_purity(self.obs, self.mc_adata, self.index, self.chunk_size)
        self.mc_adata.obs['mc_celltype'] = evaluate.get_celltype(self.obs, self.mc_adata, self.index, self.chunk_size)
        self.res.append(self.sparsity)
        self.res.append(self.mc_sparsity)
        self.res.append(self.mc_sparsity - self.sparsity)
        self.res.append(self.mc_sparsity / self.sparsity)
        self.result.loc[0] = self.res
        self.result.to_csv(self.result_dir + '/result.csv')
        print(self.result)
        self.mc_adata.write(self.result_dir + '/%s_carousel_result.h5ad' % self.data_name)

        return self.result

