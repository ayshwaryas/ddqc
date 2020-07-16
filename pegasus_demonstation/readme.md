### This directory contains all relevant files for implementation of the data-driven filtering approach in Pegasus.

Desription of the files:
- filtering_notebook.ipynb: contains all functions for data-driven aproach
    - cluster_data: function that performs clustering; does dimensional reductions and finds DE genes if specified.
    -mad: calculates median absolute deviation on a numpy array. Constant is to make it behave similar to R MAD
    - calculate_percent_ribo: calculates percent ribo for adata. This code is taken from PG calculation of percent mito
    - initial_qc: basic qc that is performed for each method
    - metric_filter: function that performs filtering on the specified metric
    - filter_cells: wrapper function that filters data for the metrics selected
    - res: parameter for clustering resolution
    - method: parameter for filtering method (none - no additional filtering; cutoff - min 200 genes, max 10% mito; outlier and mad - data driven methods)
    - param: parameter for filtering method parameter (only relevant for mad, determines how many MADs from median is allowed)

