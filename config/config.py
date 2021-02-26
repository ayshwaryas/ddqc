from local_config import local


if not local:
    DATA_DIR = "/ahg/regevdata/projects/scqc/data/"  # directory with all data
    OUTPUT_DIR = "/ahg/regevdata/projects/scqc/output_pg/"  # output directory
else:  # for debug outside of cluster
    LOC = "cluster"

    if LOC == "local":
        DATA_DIR = "/Users/michaelalperovich/Documents/primes_storage/data/"
    elif LOC == "drive":
        DATA_DIR = "/Volumes/easystore/primes_storage/data/"
    elif LOC == "cluster":
        DATA_DIR = "/Volumes/scqc/data/"
    OUTPUT_DIR = "/Users/michaelalperovich/Documents/primes_storage/output_pg/"

SOURCE_DIR_PREFIX = OUTPUT_DIR  # for fc plots, do not change

basic_genes_filter = 100  # basic nGenes filter (performed for all methods)
basic_mito_filter = 80
resolution = 1.4  # this resolution gives results closest to seurat

mito_prefixes = {"human": "MT-", "mouse": "Mt-"}  # prefixes of mitochondrial genes
# TODO: check ribo prefix
ribo_prefixes = {"human": "^Rp[sl]\d", "mouse": "^Rp[sl]\d"}  # prefixes of ribosomal genes

# if true, filtering will be done for the selected metric
do_counts = True
do_genes = True
do_mito = True
do_ribo = True
