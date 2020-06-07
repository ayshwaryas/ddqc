from local_config import local


if not local:
    DATA_DIR = "/ahg/regevdata/projects/scqc/data/"
    OUTPUT_DIR = "/ahg/regevdata/projects/scqc/output/"
else:
    DATA_DIR = "/Users/michaelalperovich/Documents/primes_storage/data/"
    OUTPUT_DIR = "/Users/michaelalperovich/Documents/primes_storage/output_pg/"
SOURCE_DIR_PREFIX = OUTPUT_DIR

CELLS_FILTER = 3
FEATURES_FILTER = 100

do_counts = True
do_genes = True
do_mito = True
do_ribo = True
