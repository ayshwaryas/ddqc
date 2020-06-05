import numpy as np
import pegasus as pg

import paths
from utils import cluster_data

INF = 10 ** 10


def mad(data, axis=0, constant=1.4826):
    return constant * np.median(np.absolute(data - np.median(data, axis)), axis)


def calculate_percent_ribo(adata, ribo_prefix):
    ribo_prefixes = ribo_prefix.split(",")

    def startswith(name):
        for prefix in ribo_prefixes:
            if name.startswith(prefix):
                return True
        return False

    ribo_genes = adata.var_names.map(startswith).values.nonzero()[0]
    adata.obs["percent_ribo"] = (adata.X[:, ribo_genes].sum(axis=1).A1 / np.maximum(adata.obs["n_counts"].values,
                                                                                    1.0)) * 100
    return adata


def initial_qc(adata, n_genes, n_cells):
    pg.qc_metrics(adata, mito_prefix="mt-", min_umis=-INF, max_umis=INF, min_genes=n_genes, max_genes=INF,
                  percent_mito=80)
    adata = calculate_percent_ribo(adata, "Rpl,Rps")
    adata = adata[:, adata.var.n_cells >= n_cells]
    pg.filter_data(adata)
    return adata


def metric_filter(adata, method, param, metric_name, do_lower_co=False, do_upper_co=False, param2=None):
    adata.obs[metric_name + "_qc_pass"] = False
    for cl in range(1, max(list(map(int, adata.obs.louvain_labels.cat.categories))) + 1):
        lower_co = -INF
        upper_co = INF
        cluster = adata[adata.obs.louvain_labels == str(cl)]
        if method == "mad":
            if do_lower_co:
                lower_co = np.median(cluster.obs[metric_name]) - param * mad(cluster.obs[metric_name])
            if do_upper_co:
                upper_co = np.median(cluster.obs[metric_name]) + param * mad(cluster.obs[metric_name])
        if method == "outlier":
            q75, q25 = np.percentile(cluster.obs[metric_name], [75, 25])
            if do_lower_co:
                lower_co = q25 - 1.5 * (q75 - q25)
            if do_upper_co:
                upper_co = q75 + 1.5 * (q75 - q25)
        if method == "cutoff":
            if do_lower_co and do_upper_co:
                lower_co = param
                upper_co = param2
            if do_lower_co:
                lower_co = param
            if do_upper_co:
                upper_co = param
        filters = [
            adata.obs.louvain_labels == str(cl),
            adata.obs[metric_name] >= lower_co,
            adata.obs[metric_name] <= upper_co
        ]
        adata.obs.loc[np.logical_and.reduce(filters), metric_name + "_qc_pass"] = True
    with open(paths.results_dir + "!filtered_" + metric_name[metric_name.find("_") + 1:] + ".csv", "w") as file:
        file.write(adata[adata.obs[metric_name + "_qc_pass"] == False].obs.to_csv())
    return adata


def filter_cells(adata, res, method, threshold, do_counts, do_genes, do_mito, do_ribo):
    adata = initial_qc(adata, 100, 3)
    adata_copy = adata[:]
    adata_copy = cluster_data(adata_copy, res)

    if method == "none":
        return adata_copy
    if do_counts:
        adata_copy = metric_filter(adata_copy, method, threshold, "n_counts", do_lower_co=True)
    else:
        adata_copy.obs["n_counts_qc_pass"] = True
    if do_genes:
        adata_copy = metric_filter(adata_copy, method, threshold, "n_genes", do_lower_co=True)
    else:
        adata_copy.obs["n_genes_qc_pass"] = True
    if do_mito:
        adata_copy = metric_filter(adata_copy, method, threshold, "percent_mito", do_upper_co=True)
    else:
        adata_copy.obs["percent_mito_qc_pass"] = True
    if do_ribo:
        adata_copy = metric_filter(adata_copy, method, threshold, "percent_ribo", do_upper_co=True)
    else:
        adata_copy.obs["percent_ribo_qc_pass"] = True
    filters = [
        adata_copy.obs["n_counts_qc_pass"],
        adata_copy.obs["n_genes_qc_pass"],
        adata_copy.obs["percent_mito_qc_pass"],
        adata_copy.obs["percent_ribo_qc_pass"],
    ]
    adata_copy.obs["passed_qc"] = False
    adata_copy.obs.loc[np.logical_and.reduce(filters), "passed_qc"] = True
    adata.obs["passed_qc"] = adata_copy.obs.passed_qc
    pg.filter_data(adata)
    return adata
