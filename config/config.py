from config.local_config import local


LOCAL = local  # needed for running outside the cluster

# PATHS
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

# CLUSTERING
resolution = 1.4  # this resolution gives results closest to seurat

# FILTERING
basic_genes_filter = 100  # basic nGenes filter (performed for all methods)
basic_mito_filter = 80
MITO_PREFIXES = {"human": "MT-", "mouse": "Mt-"}  # prefixes of mitochondrial genes
# TODO: check ribo prefix
RIBO_PREFIXES = {"human": "^Rp[sl]\d", "mouse": "^Rp[sl]\d"}  # prefixes of ribosomal genes

# if true, filtering will be done for the selected metric
do_counts = True
do_genes = True
do_mito = True
do_ribo = True

# METHOD COMPARISON SCRIPTS
# none - no additional filtering; cutoff - min 200 genes, max 10% mito; outlier and mad - data driven methods
MC_METHODS = (("none", 0), ("cutoff", 10), ("mad", 2))
MC_TASKS_PER_TISSUE = len(MC_METHODS)
