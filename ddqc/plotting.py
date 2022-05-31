from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pegasusio import UnimodalData, MultimodalData

from ddqc.filtering import perform_ddqc


def boxplot_sorted(df: pd.DataFrame, column: str, by: str, hline_x: Union[int, None] = None, log: bool = False):
    """Function that plots ddqc boxplot for specified metric."""
    df_copy = df.loc[:, [column, by]]
    if log:
        df_copy[column] = np.log2(df_copy[column])
    my_order = df_copy.groupby(by=by)[column].median().sort_values(ascending=False).index
    fig, ax = plt.subplots(figsize=(15, 12))
    chart = sns.boxplot(ax=ax, x=df_copy[by], y=df_copy[column], order=my_order)
    if log:
        ax.set_ylabel("log2({})".format(column))
    if hline_x is not None:
        chart.axhline(hline_x, color="red")
    return chart


def calculate_filtering_stats(data: MultimodalData, threshold: float, n_genes_lower_bound: int,
                              percent_mito_upper_bound: float):
    """Helper function that calculated filtering stats for filtering_facet_plot."""
    start, end = min(10, int(threshold * 10)), max(30, int(threshold * 10))
    thresholds = []
    while start <= end:
        if start != threshold:
            thresholds.append(start / 10)
        start += 1
    thresholds.append(threshold)
    filtering_stats = pd.DataFrame(columns=["threshold", "cluster", "filtered_cells", "filtered_cells_pct"])

    for th in thresholds:
        _, _, filtering_stats = perform_ddqc(data, "mad", th, 0, 0, 0, 0, n_genes_lower_bound,
                                             percent_mito_upper_bound, filtering_stats)
    return filtering_stats


def filtering_facet_plot(plot_data: pd.DataFrame, threshold: float, pct=False) -> None:
    """Function that plots facet plot of number of cells filtered out for threshold per cluster."""
    if not pct:
        y_col_name = "filtered_cells"
        y_ax_label = "Number of cells filtered"
    else:
        y_col_name = "filtered_cells_pct"
        y_ax_label = "Percentage of cells filtered"
    grid = sns.FacetGrid(plot_data, col="cluster", col_wrap=5, sharex=False, sharey=False)
    grid.map(sns.lineplot, "threshold", y_col_name)
    grid.set_axis_labels(x_var="MAD Threshold", y_var=y_ax_label)
    grid.refline(x=threshold)
