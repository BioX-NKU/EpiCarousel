import epicarousel
import sys

data_name = sys.argv[1]
carousel_resolution = int(sys.argv[2])  # 10 or 75 for reproducibility
# index = sys.argv[3]

print(data_name)
print(carousel_resolution)
# print(index)

result_base = '/home/metacell/data/metacell/EpiCarousel/Results'
epicarousel.core.setup_seed(1)
data_dir = '/home/metacell/data/metacell/source_data/source_data_shuffled_seed_1/%s_shuffled.h5ad'%data_name
if_bi = 1
if_mc_bi = 1
threshold = 0
filter_rate = 0.01
chunk_size = 10000
# carousel_resolution = 10
step = 20
threads=8
if data_name == 'AHT' or data_name == 'HFO':
    index = 'cell type'
else:
    index = 'cell_type'
    
print(index)
mc_mode = 'average'
carousel = epicarousel.core.Carousel(data_name=data_name,
                                     data_dir=data_dir,
                                     if_bi=if_bi,
                                     if_mc_bi=if_mc_bi,
                                     threshold=threshold,
                                     filter_rate=filter_rate,
                                     chunk_size=chunk_size,
                                     carousel_resolution=carousel_resolution,
                                     base=result_base,
                                     step=step,
                                     threads=threads,
                                     index=index,
                                     mc_mode=mc_mode
                                    )
carousel.make_dirs()
carousel.data_split()
carousel.identify_metacells()
carousel.merge_metacells()
carousel.metacell_preprocess()
carousel.metacell_data_clustering()
carousel.result_comparison()
carousel.delete_dirs()