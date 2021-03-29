import os
from collections import Counter

import numpy as np
import pandas as pd
import pegasus as pg

from config.config import DATA_DIR, do_batch_correction


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


# convert marker_dict to pandas df
def marker_dict_to_df(marker_dict, min_log_fc=0.25, min_pct=25, max_pval=0.05):
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
    markers = pd.concat(frames)  # merge all marker data frames together
    # sort markers by cluster, then by Log2FC
    markers.sort_values(by=["cluster", "log2FC"], ascending=[True, False], inplace=True)
    return markers


# function that performs clustering; does dimensional reductions and finds DE genes if specified
def cluster_data(adata, resolution, compute_markers=False, compute_reductions=False):
    pg.log_norm(adata)
    if do_batch_correction:
        pg.highly_variable_features(adata, consider_batch=True)
        pg.pca(adata, random_state=29)
        pca_key = pg.run_harmony(adata, random_state=29)
        pg.neighbors(adata, K=20, rep=pca_key, random_state=29)
        pg.louvain(adata, resolution=resolution, rep=pca_key, random_state=29)
    else:
        pg.highly_variable_features(adata, consider_batch=False)
        pg.pca(adata, random_state=29)
        pg.neighbors(adata, K=20, random_state=29)
        pg.louvain(adata, resolution=resolution, random_state=29)

    if compute_reductions:
        pg.umap(adata, random_state=29)

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
    signature_dir = DATA_DIR + "signatures/"
    signature_list = pd.read_csv(signature_dir + "signatures_info.csv")
    signatures = dict()
    for i, row in signature_list.iterrows():
        if is_human and row["do_human"]:
            if row["name"].startswith("pg_"):
                pg.calc_signature_score(adata, row["file"])
            else:
                signatures[row["name"]] = [t.strip().upper() for t in open(signature_dir + row["file"]).readlines()]
        elif not is_human and row["do_mouse"]:
            if row["name"].startswith("pg_"):
                pg.calc_signature_score(adata, row["file"])
            else:
                signatures[row["name"]] = [title(t.strip()) for t in open(signature_dir + row["file"]).readlines()]
    pg.calc_signature_score(adata, signatures)
    return adata


def format_markers(markers, n=25):  # take top 25 markers from the list, and separate them with ;
    if len(markers) > 0:
        return "; ".join(markers[:min(n, len(markers))])
    else:  # if markers list is empty return empty string
        return [""]


def cell_type_replacements(cell_type, tissue):
    if cell_type == "Hepatocytes" and "liver" not in tissue.lower():
        return "Epithelial cells"
    elif cell_type == "Hepatic stellate cells" and "liver" not in tissue.lower():
        return "Epithelial cells"
    elif cell_type == "Keratinocytes" and "skin" not in tissue.lower():
        return "Epithelial cells"
    elif cell_type == "Luminal epithelial cells":
        return "Epithelial/Luminal epithelial cells"
    else:
        return cell_type


# assign cell types, annotations, and calculate cluster statistics
# min_sg is how many score genes are required to call the cell type for the cluster
def assign_cell_types(adata, markers, tissue, min_sg=3):
    print("Assigning annotations and cell types to clusters")
    genes = pd.read_csv(DATA_DIR + "markers.tsv", sep="\t")  # read panglaoDB genes to cell types table
    colnames = ["cluster", "annotation", "annotation2", "%annotation", "%annotation2", "cell_type", "cell_type2",
                "n_cells", "genes_mean", "genes_median", "mito_mean", "mito_median", "ribo_mean", "ribo_median",
                "markers", "score_genes"]  # colnames for cluster info csv
    clusters_dict = dict()  # create a dict with colnames, that will be later converted to a df
    for c in colnames:
        clusters_dict[c] = []

    for cl in range(1, max(list(map(int, adata.obs.louvain_labels.cat.categories))) + 1):  # iterate though all clusters
        cluster = adata[adata.obs.louvain_labels == str(cl)]  # subset adata based on cluster number
        cluster_markers = markers[markers.cluster == str(cl)]  # subset markers for current cluster

        possible_cell_types = dict()  # dict with possible cell types and scores (sum of Log2FC for each marker)
        score_genes = dict()  # dict with possible cell types and score_genes for them
        for index, marker_info in cluster_markers.iterrows():  # iterate through markers
            marker = index.upper()
            gene_info = genes[genes["official gene symbol"] == marker]  # find marker in PanglaoDB table
            for gene, info in gene_info.iterrows():  # iterate through  cell types for the marker
                ct = cell_type_replacements(info["cell type"], tissue)
                if ct not in possible_cell_types:  # if we didn't encounter this cell type yet, add it to dicts
                    possible_cell_types[ct] = 0
                    score_genes[ct] = []
                possible_cell_types[ct] += marker_info["log2FC"]  # increase score by Log2FC of current marker
                score_genes[ct].append(marker)  # add marker to the list of score genes for current cell type

        called = False  # whether the cell type for the cluster was called
        if len(possible_cell_types) > 1:  # more than 2 cell types
            # sort cell types by score
            srt_cell_types = sorted(list(possible_cell_types.keys()), key=lambda x: -possible_cell_types[x])
            sg = score_genes[srt_cell_types[0]]
            if len(sg) >= min_sg:  # if number of score genes satisfies the minimum requirement call the cell type
                clusters_dict["cell_type"].append(srt_cell_types[0])
                clusters_dict["cell_type2"].append(srt_cell_types[1])
                clusters_dict["score_genes"].append(format_markers(sg))
                called = True
        elif len(possible_cell_types) == 1:  # only one cell type
            ct = list(possible_cell_types.keys())[0]
            sg = score_genes[ct]
            if len(sg) >= min_sg:  # if number of score genes satisfies the minimum requirement call the cell type
                clusters_dict["cell_type"].append(ct)
                clusters_dict["cell_type2"].append("")  # second common cell type doesn't exist
                clusters_dict["score_genes"].append(format_markers(sg))
                called = True
        if not called:  # if cell type wasn't called, mark it Unknown
            clusters_dict["cell_type"].append("Unknown")
            clusters_dict["cell_type2"].append("")
            clusters_dict["score_genes"].append([])

        # calculate first and second most frequent annotation and percentage of cells that have them
        cl_anno = list(cluster.obs.annotations)
        most_common = Counter(cl_anno).most_common()
        clusters_dict["annotation"].append(most_common[0][0])
        clusters_dict["%annotation"].append(round(most_common[0][1] / len(cl_anno), 3))
        if len(most_common) > 1:  # more than one annotation
            clusters_dict["annotation2"].append(most_common[1][0])
            clusters_dict["%annotation2"].append(round(most_common[1][1] / len(cl_anno), 3))
        else:  # only one annotation
            clusters_dict["annotation2"].append("")
            clusters_dict["%annotation2"].append(0)

        # record cluster statistics and markers
        clusters_dict["cluster"].append(cl)
        clusters_dict["n_cells"].append(len(cl_anno))
        clusters_dict["genes_mean"].append(round(float(np.mean(cluster.obs.n_genes)), 3))
        clusters_dict["genes_median"].append(round(float(np.median(cluster.obs.n_genes)), 3))
        clusters_dict["mito_mean"].append(round(float(np.mean(cluster.obs.percent_mito)), 3))
        clusters_dict["mito_median"].append(round(float(np.median(cluster.obs.percent_mito)), 3))
        clusters_dict["ribo_mean"].append(round(float(np.mean(cluster.obs.percent_ribo)), 3))
        clusters_dict["ribo_median"].append(round(float(np.median(cluster.obs.percent_ribo)), 3))

        clusters_dict["markers"].append(format_markers(list(cluster_markers.index)))

        print("Cluster {} finished".format(cl))

    return pd.DataFrame.from_dict(clusters_dict)  # convert dict to pandas df
