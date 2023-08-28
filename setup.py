#!/usr/bin/env python
#-*- coding:utf-8 -*-


from setuptools import setup, find_packages

setup(
    name="epicarousel",
    version="0.0.1",
    keywords=("pip", "epicarousel"),
    description="EpiCarousel: memory- and time-efficient identification of metacells for atlas-level single-cell chromatin accessibility data",
    long_description="EpiCarousel provides an memory- and time-efficient way to identify metacells of atlas-level single-cell chromatin accessibility sequencing data for downstream analysis. We provide documentation in the form of functional application programming interface documentation, tutorials and example workflows. All EpiCarousel wheels distributed on PyPI are MIT licensed.",
    license="MIT Licence",
    url="https://github.com/BioX-NKU/EpiCarousel",
    author="Shengquan Chen, Sijie Li and Yuxi Li",
    packages=find_packages(),
    python_requires='>3.8.0',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    install_requires=[
        "scanpy>=1.8.2",
        "anndata==0.8.0",
        "snapatac2>=2.1.0",
        "episcanpy==0.3.2",
        "seaborn",
        "igraph",
        "leidenalg",
        "louvain",
        'shutilwhich',
        'memory_profiler',
        'scipy',
        'tqdm'
    ],
    package_data={
        'epicarousel': ['construction.sh']
    }
)
    
    