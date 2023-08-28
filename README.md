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

* **file path**: An h5ad file where the scCAS data count matrix is stored in the compressed sparse row format. AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to peaks/regions.
* **chunk size**:  The number of cells in each chunk.
* **resolution**: The ratio of the number of cells to that of metacells.
* **threads**:  The number of parallel processes.

