from local_config import local


if not local:
    DATA_DIR = "/ahg/regevdata/projects/scqc/data/"  # directory with all data
    OUTPUT_DIR = "/ahg/regevdata/projects/scqc/output_pg/"  # output directory
else:  # for debug outside of cluster
    DATA_DIR = "/Users/michaelalperovich/Documents/primes_storage/data/"
    DATA_DIR = "/Volumes/easystore/primes_storage/data/"
    DATA_DIR = "/Volumes/scqc/data/"
    OUTPUT_DIR = "/Users/michaelalperovich/Documents/primes_storage/output_pg/"
SOURCE_DIR_PREFIX = OUTPUT_DIR  # for fc plots, do not change

CELLS_FILTER = 3  # basic nCells filter (performed for all methods)
FEATURES_FILTER = 100  # basic nGenes filter (performed for all methods)

# if true, filtering will be done for the selected metric
do_counts = True
do_genes = True
do_mito = True
do_ribo = True
