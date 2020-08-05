import re

import numpy as np
import pegasus as pg

import paths
from utils import cluster_data

INF = 10 ** 10  # infinity for cases no filtering is required


# MAD for numpy array. Constant is to make it behave similar to seurat MAD
def mad(data, axis=0, constant=1.4826):
    return constant * np.median(np.absolute(data - np.median(data, axis)), axis)


# calculates percent ribo for adata. This code is taken from PG calculation of percent mito
def calculate_percent_ribo(adata, ribo_prefix):
    ribo_prefixes = ribo_prefix.split(",")  # parse ribo prefixes

    def startswith(name):
        for prefix in ribo_prefixes:
            if re.match(prefix, name, flags=re.IGNORECASE):
                return True
        return False

    ribo_genes = adata.var_names.map(startswith).values.nonzero()[0]  # get all genes that match the pattern
    # calculate percent ribo
    adata.obs["percent_ribo"] = (adata.X[:, ribo_genes].sum(axis=1).A1 / np.maximum(adata.obs["n_counts"].values,
                                                                                    1.0)) * 100
    return adata


# basic qc that is performed for each method
def initial_qc(adata, n_genes, n_cells, is_human):
    mt_prefix = "MT-" if is_human else "mt-"
    pg.qc_metrics(adata, mito_prefix=mt_prefix, min_umis=-INF, max_umis=INF, min_genes=n_genes, max_genes=INF,
                  percent_mito=80)  # default PG filtering with custom cutoffs
    adata = calculate_percent_ribo(adata, "^Rp[sl]\d")  # calculate percent ribo
    adata = adata[:, adata.var.n_cells > n_cells]  # filtering based on nCells
    pg.filter_data(adata)  # filtering based on the parameters from qc_metrics
    return adata


# function that performs filtering on the specified metric
# method - method name for filtering (mad, outlier, cutoff)
# param - parameter for the selected method
# metric name - name of the metric (must be in adata.obs)
# do_upper_co and do_lower_co - whether to do upper and lower cutoff
def metric_filter(adata, method, param, metric_name, do_lower_co=False, do_upper_co=False):
    adata.obs[metric_name + "_qc_pass"] = False  # T/F array to tell whether the cell is filtered
    for cl in range(1, max(list(map(int, adata.obs.louvain_labels.cat.categories))) + 1):  # iterate though all clusters
        lower_co = -INF
        upper_co = INF
        cluster = adata[adata.obs.louvain_labels == str(cl)]  # subset adata based on cluster number
        if method == "mad":  # calculate MAD cutoffs, which are median ± param * MAD
            if do_lower_co:
                lower_co = np.median(cluster.obs[metric_name]) - param * mad(cluster.obs[metric_name])
            if do_upper_co:
                upper_co = np.median(cluster.obs[metric_name]) + param * mad(cluster.obs[metric_name])
        if method == "outlier":  # calculate Outlier cutoffs, which are Q1 - 1.5 * IQR or Q3 + 1.5 * IQR
            q75, q25 = np.percentile(cluster.obs[metric_name], [75, 25])
            if do_lower_co:
                lower_co = q25 - 1.5 * (q75 - q25)
            if do_upper_co:
                upper_co = q75 + 1.5 * (q75 - q25)
        if method == "cutoff":  # cutoff uses the following filtering criteria: >= 200 genes, < 10% mito
            if metric_name == "n_counts":
                lower_co = -INF
                upper_co = INF
            if metric_name == "n_genes":
                lower_co = 200
                upper_co = INF
            if metric_name == "percent_mito":
                upper_co = 10.0
        filters = [
            adata.obs.louvain_labels == str(cl),
            adata.obs[metric_name] >= lower_co,
            adata.obs[metric_name] <= upper_co
        ]  # filtering condition
        # for cells that satisfy the condition set the value to true
        adata.obs.loc[np.logical_and.reduce(filters), metric_name + "_qc_pass"] = True
    with open(paths.results_dir + "!filtered_" + metric_name[metric_name.find("_") + 1:] + ".csv", "w") as file:
        # write the cells that failed the filtering to the csv
        file.write(adata[adata.obs[metric_name + "_qc_pass"] == False].obs.to_csv())
    return adata


# function that filters data for the metrics selected
# method - method name for filtering (mad, outlier, cutoff)
# threshold - parameter for the selected method
# do_metric - set to true, if you want to filter the data based on metric
def filter_cells(adata, res, method, threshold, is_human, do_counts, do_genes, do_mito, do_ribo):
    adata = initial_qc(adata, 100, 3, is_human)  # perform initial qc with min 100 genes and min 3 cells
    adata_copy = adata.copy()  # make a copy of adata, so the clustering results won't affect future downstream analysis
    adata_copy = cluster_data(adata_copy, res)  # do initial clustering of adata

    if method == "none":  # no filtering option
        return adata

    # for each metric if do_metric is true, the filtering will be performed
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
    ]  # filtering condition
    adata_copy.obs["passed_qc"] = False  # cumulative T/F array for to determine cells that passed the filtering
    # for cells that satisfy the condition set the value to true
    adata_copy.obs.loc[np.logical_and.reduce(filters), "passed_qc"] = True
    adata.obs["passed_qc"] = adata_copy.obs.passed_qc  # transfer array from the copy to actual object
    pg.filter_data(adata)  # perform filtering
    return adata
