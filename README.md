# EpiCarousel: memory- and time-efficient identification of metacells for atlas-level single-cell chromatin accessibility data



## Installation
EpiCarousel is available on PyPI and can be installed via

```
pip install epicarousel
```

You can also install epicarousel from GitHub via
```
git clone git://github.com/BioX-NKU/EpiCarousel.git
cd EpiCarousel
python setup.py install
```
The dependencies will be automatically installed along with EpiCarousel.


## Quick Start

### Input

* **data_name**: Name of the dataset.
* **data_dir**: Path of an h5ad file where the scCAS data count matrix is stored in the compressed sparse row format. AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to peaks/regions.
* **if_bi**: Whether to binarize the scCAS data count matrix.
* **if_mc_bi**: Whether to binarize the metacell-by-region matrix. 
* **threshold**: Threshold for binarizing metacell-by-region matrix.
* **filter_rate**: Proportion for feature selection.
* **chunk_size**:  Number of cells in each chunk.
* **carousel_resolution**: Ratio of the number of cells to that of metacells.
* **base**: Export path for EpiCarousel.
* **step**: Length of Walktrap community detection.
* **threads**: Number of parallel processes.
* **mc_mode**: Mode of calculating metacell-by-region matrix.
* **index**: (Optional) Ground truth cell type label of single cells for downstream analysis and evaluation.

### Output

+ adata: Metacell AnnData object of shape `n_obs` × `n_vars` stored in an h5ad file. Rows correspond to metacells and columns to features.

EpiCarousel can also be seamlessly integrated with [EpiScanpy](https://episcanpy.readthedocs.io/en/stable/), a widely-used Python library for epigenomics single cell analysis:

```Python
import episcanpy.api as epi
import epicarousel

# Run EpiCarousel
carousel = epicarousel.core.Carousel(data_name,
                                     data_dir,
                                     if_bi,
                                     if_mc_bi,
                                     threshold,
                                     filter_rate,
                                     chunk_size,
                                     carousel_resolution,
                                     base,
                                     step,
                                     threads,
                                     mc_mode,
                                     index
                                    )
carousel.make_dirs()
carousel.data_split()
carousel.metacell_construction()
carousel.merge_metacell()
carousel.metacell_preprocess()
carousel.metacell_data_clustering()
carousel.result_comparison()
carousel.delete_dirs()

# Load the metacell data as an AnnData object (adata).
```

## The source code for the reproduction of results can be found [here](https://github.com/BioX-NKU/EpiCarousel_reporducibility/).

## We also provide a [Jupyter Notebook]() for running EpiCarousel.

## Find more details on [the Documentation of EpiCarousel]().





