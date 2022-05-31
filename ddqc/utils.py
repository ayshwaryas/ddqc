from typing import Tuple

import numpy as np
import pandas as pd
import pegasus as pg
from pegasusio import UnimodalData, MultimodalData


def mad(x: np.ndarray, constant: float = 1.4826) -> float:
    """Function that computes adjusted MAD for numpy array"""
    return constant * np.median(np.absolute(x - np.median(x)))


def calculate_percent_ribo(data: MultimodalData, ribo_prefix) -> None:
    """Function that calculates percent ribo for Pegasus object"""
    import re
    ribo_genes = \
        data.var_names.map(lambda x: re.match(ribo_prefix, x, flags=re.IGNORECASE) is not None).values.nonzero()[
            0]  # get all genes that match the pattern
    data.obs["percent_ribo"] = (data.X[:, ribo_genes].sum(axis=1).A1 / np.maximum(data.obs["n_counts"].values,
                                                                                  1.0)) * 100  # calculate percent ribo


def cluster_data(data: MultimodalData, n_genes: int, percent_mito: float, mito_prefix: str, ribo_prefix: str,
                 norm_count: float = 1e5, n_components: int = 50, k: int = 20, clustering_method: str = "louvain",
                 resolution: float = 1.3, random_state: int = 29) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Function that clusters a Pegasus object before filtering"""
    pg.qc_metrics(data, mito_prefix=mito_prefix, min_umis=-10 ** 10, max_umis=10 ** 10, min_genes=n_genes,
                  max_genes=10 ** 10,
                  percent_mito=percent_mito)  # default PG filtering with custom cutoffs
    calculate_percent_ribo(data, ribo_prefix)  # calculate percent ribo
    pg.filter_data(data)  # filtering based on the parameters from qc_metrics
    pg.identify_robust_genes(data)

    obs_copy = data.obs.copy()
    var_copy = data.var.copy()
    var_copy.drop(columns=["robust", "highly_variable_features", "n_cells", "percent_cells"], inplace=True,
                  errors="ignore")
    uns_copy = data.uns.mapping.copy()

    pg.log_norm(data, norm_count=norm_count)
    pg.highly_variable_features(data, consider_batch=False)
    pg.pca(data, n_components=n_components, random_state=random_state)
    pg.neighbors(data, K=k, random_state=random_state)
    if clustering_method == "louvain":
        pg.louvain(data, resolution=resolution, random_state=random_state)
    elif clustering_method == "leiden":
        pg.leiden(data, resolution=resolution, random_state=random_state)
    elif clustering_method == "spectral_louvain":
        pg.spectral_louvain(data, resolution=resolution, random_state=random_state)
    elif clustering_method == "spectral_leiden":
        pg.spectral_leiden(data, resolution=resolution, random_state=random_state)
    else:
        raise KeyError(f"Unknown clustering method {clustering_method}")

    data.obs["cluster_labels"] = data.obs[clustering_method + "_labels"]
    data.obs[clustering_method + "_labels"] = None
    return obs_copy, var_copy, uns_copy


def reverse_to_raw_matrix(unidata: UnimodalData, obs_copy: pd.DataFrame, var_copy: pd.DataFrame, uns_copy: dict):
    """Function that reverses a Pegasus object to a raw matrix and removes all additional information"""
    unidata.obs = obs_copy
    unidata.var = var_copy
    unidata.matrices["X"] = unidata.matrices.pop("raw.X")
    unidata.obsm.clear()
    unidata.varm.clear()
    unidata.uns = uns_copy
