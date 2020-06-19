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


# currently not used
def load_markers_dict():
    markers_dict = dict()
    markers = pd.read_csv(DATA_DIR + "markers.tsv", delimiter="\t")

    for index, row in markers.iterrows():
        markers_dict[row["official gene symbol"]] = row["cell type"]
    return markers_dict
