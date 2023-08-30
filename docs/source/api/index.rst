.. module:: epicarousel
.. automodule:: epicarousel
   :noindex:
API
====


Import EpiCarousel::

   import epicarousel
   
Carousel
------------
.. module::epicarousel.core
.. currentmodule::epicarousel

.. autosummary::
    :toctree: .

    core.Carousel

Data processing
-----------
.. module::epicarousel.core
.. currentmodule::epicarousel


.. autosummary::
    :toctree: .

    core.Carousel.make_dirs
    core.Carousel.data_split
    core.Carousel.delete_dirs
    
Metacell identification
----------
.. module::epicarousel.core
.. currentmodule::epicarousel

.. autosummary::
    :toctree: .
    
    core.Carousel.identify_metacells
    
Metacell aggregation
-----------
.. module::epicarousel.core
.. currentmodule::epicarousel


.. autosummary::
    :toctree: .

    core.Carousel.merge_metacells

Downstream analysis
-----------
.. module::epicarousel.core
.. currentmodule::epicarousel


.. autosummary::
    :toctree: .

    core.Carousel.metacell_preprocess
    core.Carousel.metacell_data_clustering

Evaluation
-----------
.. module::epicarousel.core
.. currentmodule::epicarousel


.. autosummary::
    :toctree: .

    core.Carousel.result_comparison
