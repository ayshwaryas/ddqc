import os

import pandas as pd
import pegasus as pg

from config.config import DATA_DIR


# check if the dir exists, if not - create it
def safe_mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)


# save all relevant cell info to csv for plots in seurat
def save_to_csv(adata, path):
    df = adata.obs

    # add dimensional reductions to data frame
    df["pca1"] = [t[0] for t in list(adata.obsm["X_pca"])]
    df["pca2"] = [t[1] for t in list(adata.obsm["X_pca"])]

    # tsne is turned off for now
    df["umap1"] = [t[0] for t in list(adata.obsm["X_umap"])]
    df["umap2"] = [t[1] for t in list(adata.obsm["X_umap"])]

    with open(path + "!cells.csv", "w") as fout:
        fout.write(df.to_csv())  # write df to csv


# write markers to csv
def write_markers(marker_dict, path, min_log_fc=0.25, min_pct=25, max_pval=0.05):
    frames = []
    # iterate through all keys in the markers dict
    for cl in marker_dict.keys():
        for d in marker_dict[cl].keys():
            df = marker_dict[cl][d]
            df['cluster'] = cl
            df['up/down'] = d

            # filter markers based on log_fc and pct
            df = df[(df["t_pval"] <= max_pval) & (df["log2FC"] >= min_log_fc) & (
                    (df["percentage"] >= min_pct) | (df["percentage_other"] >= min_pct))]

            frames.append(df)
    result = pd.concat(frames)  # merge all marker data frames together
    with open(path + "!markers.csv", "w") as fout:
        fout.write(result.to_csv())


# function that performs clustering; does dimensional reductions and finds DE genes if specified
def cluster_data(adata, resolution, compute_markers=False, compute_reductions=False):
    pg.log_norm(adata)
    pg.highly_variable_features(adata, consider_batch=False)
    pg.pca(adata)
    pg.neighbors(adata, K=20)  # TODO: should we keep K changed
    pg.louvain(adata, resolution=resolution)

    if compute_reductions:
        pg.umap(adata)

    if compute_markers:
        pg.de_analysis(adata, cluster='louvain_labels', t=True, fisher=False, temp_folder="/tmp")
        marker_dict = pg.markers(adata, alpha=1)
        return adata, marker_dict
    else:
        return adata


def title(s):  # convert string to title case
    if len(s) == 1:
        return s[0].upper()
    return s[0].upper() + s[1:].lower()


def add_cd_scores(adata, is_human):  # add cell death signatures
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
