import os

import pegasus as pg
import pandas as pd

from config import DATA_DIR


# check if the dir exists, if not - create it
def safe_mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)


# function that performs clustering; does dimensional reductions and finds DE genes if specified
def cluster_data(adata, resolution=1, compute_markers=False, compute_reductions=False):
    pg.log_norm(adata)
    pg.highly_variable_features(adata, consider_batch=False)
    pg.pca(adata)
    pg.neighbors(adata, K=20)  # K=20 to make it closer to seurat

    if resolution == 0:  # if resolution == 0, run louvain with default resolution (1.3)
        pg.louvain(adata)
    else:
        pg.louvain(adata, resolution=resolution)

    if compute_reductions:
        # pg.fitsne(adata)
        pg.umap(adata)

    if compute_markers:
        pg.de_analysis(adata, cluster='louvain_labels', auc=False, t=True, fisher=False, mwu=False,
                       temp_folder="/tmp")
        marker_dict = pg.markers(adata)
        return adata, marker_dict
    else:
        return adata


def title(s):
    if len(s) == 1:
        return s[0].upper()
    return s[0].upper() + s[1:].lower()


def add_cd_scores(adata, is_human):
    if is_human:
        cd1 = [t.strip().upper() for t in open(DATA_DIR + "signatures/cd1_signatures.csv").readlines()]
        cd2 = [t.strip().upper() for t in open(DATA_DIR + "signatures/cd2_signatures.csv").readlines()]
        cd3 = [t.strip().upper() for t in open(DATA_DIR + "signatures/cd3_signatures.csv").readlines()]
        mito = "mitochondrial_genes_human"
        ribo = "ribosomal_genes_human"
    else:
        cd1 = [title(t.strip()) for t in open(DATA_DIR + "signatures/cd1_signatures.csv").readlines()]
        cd2 = [title(t.strip()) for t in open(DATA_DIR + "signatures/cd2_signatures.csv").readlines()]
        cd3 = [title(t.strip()) for t in open(DATA_DIR + "signatures/cd3_signatures.csv").readlines()]
        mito = "mitochondrial_genes_mouse"
        ribo = "ribosomal_genes_mouse"
    signatures = {"cd1": cd1, "cd2": cd2, "cd3": cd3}

    pg.calc_signature_score(adata, signatures)
    pg.calc_signature_score(adata, mito)
    pg.calc_signature_score(adata, ribo)
    return adata


# currently not used
def load_markers_dict():
    markers_dict = dict()
    markers = pd.read_csv(DATA_DIR + "markers.tsv", delimiter="\t")

    for index, row in markers.iterrows():
        markers_dict[row["official gene symbol"]] = row["cell type"]
    return markers_dict
