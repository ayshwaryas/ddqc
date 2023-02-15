from typing import Union

import numpy as np
import pandas as pd
from pegasusio import MultimodalData

from ddqc.utils import mad


def metric_filter(data: MultimodalData, method: str, param: float, metric_name: str, do_lower_co: bool = False,
                  do_upper_co: bool = False,
                  lower_bound: float = 10 ** 10, upper_bound: float = -10 ** 10,
                  df_qc: pd.DataFrame = None) -> np.ndarray:
    """Function that performs filtering result computation for the specified metric."""
    qc_pass = np.zeros(data.shape[0], dtype=bool)  # T/F array to tell whether the cell is filtered
    if df_qc is not None:
        df_qc[f"{metric_name}_lower_co"] = None
        df_qc[f"{metric_name}_upper_co"] = None

    for cl in data.obs["cluster_labels"].cat.categories:  # iterate through all clusters
        idx = data.obs["cluster_labels"] == cl
        values = data.obs.loc[idx, metric_name]

        if method == "mad":  # calculate MAD cutoffs, which are median ± param * MAD
            median_v = np.median(values)
            mad_v = mad(values)
            lower_co = median_v - param * mad_v
            upper_co = median_v + param * mad_v
        else:
            assert method == "outlier"  # calculate Outlier cutoffs, which are Q1 - 1.5 * IQR or Q3 + 1.5 * IQR
            q75, q25 = np.percentile(values, [75, 25])
            lower_co = q25 - 1.5 * (q75 - q25)
            upper_co = q75 + 1.5 * (q75 - q25)

        lower_co = min(lower_co, lower_bound)
        upper_co = max(upper_co, upper_bound)

        qc_pass_cl = np.ones(values.size, dtype=bool)
        if df_qc is not None:
            df_qc.loc[idx, f"{metric_name}"] = values
        if do_lower_co:
            qc_pass_cl &= (values >= lower_co)
            if df_qc is not None:
                df_qc.loc[idx, f"{metric_name}_lower_co"] = lower_co
        if do_upper_co:
            qc_pass_cl &= (values <= upper_co)
            if df_qc is not None:
                df_qc.loc[idx, f"{metric_name}_upper_co"] = upper_co
        if df_qc is not None:
            df_qc.loc[idx, f"{metric_name}_passed_qc"] = qc_pass_cl
        qc_pass[idx] = qc_pass_cl
    return qc_pass


def perform_ddqc(data: MultimodalData, method: str, threshold: float,
                 threshold_counts: Union[int, None], threshold_genes: Union[int, None],
                 threshold_mito: Union[float, None], threshold_ribo: Union[float, None],
                 n_genes_lower_bound: int, percent_mito_upper_bound: float,
                 filtering_stats: Union[pd.DataFrame, None] = None):
    """Function that computes ddqc for all requested metrics with specified parameters."""
    df_qc = pd.DataFrame({"cluster_labels": data.obs["cluster_labels"].values}, index=data.obs_names)
    passed_qc = np.ones(data.shape[0], dtype=bool)

    threshold_counts = threshold if threshold_counts == 0 else threshold_counts
    threshold_genes = threshold if threshold_genes == 0 else threshold_genes
    threshold_mito = threshold if threshold_mito == 0 else threshold_mito
    threshold_ribo = threshold if threshold_ribo == 0 else threshold_ribo

    # for each metric in case do_metric is true, the filtering will be performed
    if threshold_counts is not None:
        passed_qc &= metric_filter(data, method, threshold_counts, "n_counts", do_lower_co=True, df_qc=df_qc)
    if threshold_genes is not None:
        passed_qc &= metric_filter(data, method, threshold_genes, "n_genes", do_lower_co=True,
                                   lower_bound=n_genes_lower_bound, df_qc=df_qc)
    if threshold_mito is not None:
        passed_qc &= metric_filter(data, method, threshold_mito, "percent_mito", do_upper_co=True,
                                   upper_bound=percent_mito_upper_bound, df_qc=df_qc)
    if threshold_ribo is not None:
        passed_qc &= metric_filter(data, method, threshold_ribo, "percent_ribo", do_upper_co=True, df_qc=df_qc)

    if filtering_stats is not None:
        cl_filtering_stats = []
        for cl in data.obs.cluster_labels.cat.categories:
            idx = data.obs["cluster_labels"] == cl
            unique, counts = np.unique(passed_qc[idx], return_counts=True)
            n_filtered_cells = dict(zip(unique, counts)).get(False, 0)
            n_filtered_cells_pct = n_filtered_cells / len(passed_qc[idx]) * 100
            row = {"threshold": threshold, "cluster": cl, "filtered_cells": n_filtered_cells,
                   "filtered_cells_pct": n_filtered_cells_pct}
            cl_filtering_stats.append(row)
        filtering_stats = pd.concat([filtering_stats, pd.DataFrame.from_records(cl_filtering_stats)])
    return passed_qc, df_qc, filtering_stats
