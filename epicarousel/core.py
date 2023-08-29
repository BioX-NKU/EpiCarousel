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
                 carousel_resolution: int,
                 step: int,
                 threads: int,
                 if_bi: int,
                 if_mc_bi: int,
                 mc_mode: str,
                 threshold: float,
                 filter_rate: float,
                 neighbors_method: str,
                 n_components: int,
                 svd_solver: str
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
    carousel_resolution : int
        The ratio of the number of original cells to that of metacells.
    step : int
        Walktrap length.
    threads : int
        Number of the processes.
    if_bi : int
        Whether to binarize the scCAS data count matrix.
    if_mc_bi : int
        Whether to binarize the metacell-by-region matrix. 
    mc_mode : str
        Mode of calculating metacell-by-region matrix.
    threshold : float
        Threshold for binarizing metacell-by-region matrix.
    filter_rate : float
        Proportion for feature selection.
    """        
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
              + ' ' + str(carousel_resolution)
              + ' ' + str(step)
              + ' ' + str(threads)
              + ' ' + str(if_bi)
              + ' ' + str(if_mc_bi)
              + ' ' + mc_mode
              + ' ' + str(threshold)
              + ' ' + str(filter_rate)
              + ' ' + neighbors_method
              + ' ' + str(n_components)
              + ' ' + svd_solver
             )

class Carousel:

    def __init__(self,
                 data_name: str,
                 data_dir: str,
                 if_bi: int = 1,
                 if_mc_bi: int = 1,
                 threshold: float = 0.0,
                 filter_rate: float = 0.01,
                 chunk_size: int = 10000,
                 carousel_resolution: int = 10,
                 base: str = '/home/metacell/data/metacell/carousel/output',
                 step: int = 4,
                 threads: int = 8,
                 mc_mode: str = 'average',
                 index: str = 'cell_type',
                 neighbors_method: str = 'umap',
                 n_components: int = 50,
                 svd_solver: str = 'arpack',
                 shuffle: int = 0,
                 random_state: int = 1,
                 ):
        '''

        Args:
            data_name: Name of data.
            data_dir: Path where the data is saved.
            if_bi: int, whether the original data is binarized, 1 represents binarization
                and other numbers represent no binarization
            if_mc_bi: int, whether metacell count matrix is binared, 1 represents binarization
                        and other numbers represent no binarization
            threshold: float, the threshold for metacell count matrix
            filter_rate: min_cells=np.ceil(filter_rate*adata.shape[0])
            chunk_size: Number of cells in each chunk.
            carousel_resolution: Ratio of the number of cells to that of metacells.
            base: Export path for EpiCarousel.
            step: Length of Walktrap community detection.
            threads: Number of parallel processes.
            index: (Optional) Ground truth cell type label of single cells for downstream analysis and evaluation.
            mc_mode: Mode of calculating metacell-by-region matrix.
        '''
        self.data_name = data_name
        self.chunk_size = chunk_size
        self.carousel_resolution = carousel_resolution
        self.gamma = 1 / carousel_resolution
        self.if_bi = if_bi  
        self.if_mc_bi = if_mc_bi  
        self.threshold = threshold
        self.random_state = random_state
        self.step = step
        self.threads = threads  
        self.shuffle = shuffle
        self.index = index
        self.mc_mode = mc_mode
        self.filter_rate = filter_rate
        self.neighbors_method = neighbors_method
        self.n_components = n_components
        self.svd_solver = svd_solver

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
        # Create output directories.

        '''
        self.proj_dir = self.base + '/EpiCarousel_export'
        print("Project output directory: %s" % self.proj_dir)
        self.sf_dir = self.proj_dir + '/shuffled'
        self.chunk_dir = self.proj_dir + '/chunk'
        self.chunk_preprocessed_dir = self.proj_dir + '/chunk_preprocessed'
        self.chunk_mc_dir = self.proj_dir + '/chunk_mc'
        self.chunk_mc_binarized_dir = self.proj_dir + '/chunk_mc_binarized'
        self.carousel_dir = self.proj_dir + '/carousel'
        self.carousel_binarized_dir = self.proj_dir + '/carousel_binarized'
        self.carousel_preprocessed_dir = self.proj_dir + '/carousel_preprocessed'
        self.carousel_clustered_dir = self.proj_dir + '/carousel_clustered'
        self.result_dir = self.proj_dir + '/carousel_clustering_results'

        self.dirs = [self.proj_dir, self.sf_dir, self.chunk_dir, self.chunk_preprocessed_dir, self.chunk_mc_dir,
                     self.chunk_mc_binarized_dir, self.carousel_dir, self.carousel_binarized_dir,
                     self.carousel_preprocessed_dir, self.carousel_clustered_dir, self.result_dir]
        for directory in self.dirs:
            if not os.path.isdir(directory):
                os.makedirs(directory)

#     @profile
    def shuffle_data(self):
        """
        Shuffle the initial data.
        """
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
        """
        Partition data sequentially into chunks.
        """
        print("Start spliting data!")
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
        print("Finish spliting data.")

#     @profile
    def identify_metacells(self):
        """
        Identify metacells.
        """
        constrcution(fold_number=self.fold_number,
                     data_name=self.data_name,
                     chunk_dir=self.chunk_dir,
                     chunk_preprocessed_dir=self.chunk_preprocessed_dir,
                     chunk_mc_dir=self.chunk_mc_dir,
                     chunk_mc_binarized_dir=self.chunk_mc_binarized_dir,
                     carousel_resolution=self.carousel_resolution,
                     step=self.step,
                     threads=self.threads,
                     if_bi=self.if_bi,
                     if_mc_bi=self.if_mc_bi,
                     mc_mode=self.mc_mode,
                     threshold=self.threshold,
                     filter_rate=self.filter_rate,
                     neighbors_method=self.neighbors_method,
                     n_components=self.n_components,
                     svd_solver=self.svd_solver
                     )

#     @profile
    def merge_metacells(self):
        '''
        Aggregate metacells from each chunk.
        '''

        print("Start aggregating metacells.")
        if self.if_mc_bi == 1:
            i = 0
            self.mc_adata = pp.read(self.chunk_mc_binarized_dir + "/%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (
            self.data_name, i + 1, self.if_bi, self.if_mc_bi))
            for i in tqdm(range(1, self.fold_number)):
                mc_adata = pp.read(self.chunk_mc_binarized_dir + "/%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (
                self.data_name, i + 1, self.if_bi, self.if_mc_bi))
                self.mc_adata = ad.concat([self.mc_adata, mc_adata])
                del mc_adata
                gc.collect()

            self.mc_adata.obs_names_make_unique()
            print("Finish aggregating metacells.")
            try:
                for varname in self.var[:].columns:
                    self.mc_adata.var[varname] = pd.Categorical(self.var[varname])
                del self.mc_adata.var['_index']
            except:
                pass

            # Save binarized metacell-by-region matrix.
            self.mc_adata.write(self.carousel_binarized_dir + '/%s_carousel_binarized.h5ad' % self.data_name)
        else:
            i = 0
            self.mc_adata = pp.read(self.chunk_mc_dir + "/%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (
            self.data_name, i + 1, self.if_bi, self.if_mc_bi))
            for i in tqdm(range(1, self.fold_number)):
                mc_adata = pp.read(self.chunk_mc_dir + "/%s_mc%d_if_bi_%d_if_mc_bi_%d.h5ad" % (
                self.data_name, i + 1, self.if_bi, self.if_mc_bi))
                self.mc_adata = ad.concat([self.mc_adata, mc_adata])
                del mc_adata
                gc.collect()

            self.mc_adata.obs_names_make_unique()
            print("Finish aggregating metacells.")
            try:
                for varname in self.var[:].columns:
                    self.mc_adata.var[varname] = pd.Categorical(self.var[varname])
                del self.mc_adata.var['_index']
            except:
                pass

            # Save unbinarized metacell-by-region matrix.
            self.mc_adata.write(self.carousel_dir + '/%s_carousel.h5ad' % self.data_name)


#     @profile
    def metacell_preprocess(self,
                           neighbors_method: str = 'umap'
                           ):
        """
        Preprocess metacells.
        """        

        print("Start metacell preprocessing.")

        self.mc_sparsity = np.sum(self.mc_adata.X) / self.mc_adata.X.shape[0] / self.mc_adata.X.shape[1]

        self.mc_adata = pp.ATAC_preprocess(self.mc_adata,
                                           filter_rate=0.01 * (self.mc_sparsity / self.sparsity),
                                           method=neighbors_method,
                                           n_components=self.n_components,
                                           svd_solver=self.svd_solver
                                           )

#         self.mc_adata.write(self.carousel_preprocessed_dir + '/%s_carousel_pped.h5ad' % self.data_name)
        print(self.mc_adata)
        print("Finish metacell preprocessing.")

#     @profile
    def metacell_data_clustering(self):
        """
        Cluster metacells.
        """        

        print("Start metacell clustering.")
        print(self.obs)
        epi.tl.leiden(self.mc_adata, key_added='Dleiden')
        epi.tl.louvain(self.mc_adata, key_added='Dlouvain')

        epi.tl.getNClusters(self.mc_adata, n_cluster=len(self.obs[self.index].unique()), method='leiden')
        self.mc_adata.obs['Cleiden'] = self.mc_adata.obs['leiden'].copy()

        epi.tl.getNClusters(self.mc_adata, n_cluster=len(self.obs[self.index].unique()), method='louvain')
        self.mc_adata.obs['Clouvain'] = self.mc_adata.obs['louvain'].copy()
        print("Finish metacell clustering.")
#         self.mc_adata.write(self.carousel_clustered_dir + '/%s_carousel_clustered.h5ad' % self.data_name)

#     @profile
    def delete_dirs(self):
        """
        Remove intermediate files.
        """        
        for direc in [self.sf_dir, self.chunk_dir, self.chunk_preprocessed_dir, self.chunk_mc_dir,
                      self.chunk_mc_binarized_dir,self.carousel_clustered_dir]:
            rmtree(direc)

#     @profile
    def result_comparison(self):
        """
        Evaluation using four clustering strategies.
        """
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

        cluster_index = ['AMI_Dlouvain', 'ARI_Dlouvain', 'NMI_Dlouvain', 'Homo_Dlouvain', 'CS_Dlouvain', 'Vms_Dlouvain',
                         'FMI_Dlouvain',
                         'AMI_Dleiden', 'ARI_Dleiden', 'NMI_Dleiden', 'Homo_Dleiden', 'CS_Dleiden', 'Vms_Dleiden',
                         'FMI_Dleiden',
                         'AMI_Clouvain', 'ARI_Clouvain', 'NMI_Clouvain', 'Homo_Clouvain', 'CS_Clouvain', 'Vms_Clouvain',
                         'FMI_Clouvain',
                         'AMI_Cleiden', 'ARI_Cleiden', 'NMI_Cleiden', 'Homo_Cleiden', 'CS_Cleiden', 'Vms_Cleiden',
                         'FMI_Cleiden',
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

