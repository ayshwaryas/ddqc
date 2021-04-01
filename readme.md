### filtering.py contains all code for Data-Driven QC

####Description of functions:
 - cluster_data: performs default clustering(log normalization, HVG, PCA, neighbours, and louvain). Parameters:
   - adata: AnnData or MultimodalData Object
   - resolution: clustering resolution used in louvain
 - mad: finds Median Absolute Deviation for a numpy array
 - calculate_percent_ribo: calculates percent_ribo similar to how pegasus calculates percent mito. Parameters:
   - adata: AnnData or MultimodalData Object
   - ribo_prefix: string of regular expressions for ribosomal genes, separated by comma
 - initial_qc: does basic threshold QC based on number of genes and percent_mito. Also, identifies robust genes. Parameters:
   - adata: AnnData or MultimodalData Object
   - n_genes: lower threshold for n_genes
   - percent_mito: upper threshold for percent_mito
   - mito_prefix: mitochondrial genes prefix
   - ribo_prefix: ribosomal genes prefix
 - metric_filter: performs filtering on a specified metric. Data must be clustered and metric must exist in obs. Parameters:
   - adata: AnnData or MultimodalData Object
   - method: method name for filtering (mad, outlier, cutoff)
   - param: parameter for the selected method
   - metric_name: name of the metric (must be in adata.obs)
   - do_upper_co and do_lower_co: whether to do upper and lower cutoff (default: False)
   - record_path: path for recording filtered cells CSVs (keep it None if not needed) (default: None)
 - filter_cells: performs DDQC for selected metrics. Parameters:
   - adata: AnnData or MultimodalData Object
   - res: clustering resolution (default 1.3)  
   - method: method name for filtering (mad, outlier, cutoff) (default mad)
   - threshold: parameter for the selected method (default 2)
   - basic_n_genes: cutoff for basic filtering of n_genes (default 100)
   - basic_percent_mito: cutoff for basic filtering of n_genes (default 80)
   - mito_prefix: mitochondrial genes prefix (default "MT-")
   - ribo_prefix: ribosomal genes prefix (default "^Rp[sl]\d")
   - do_metric - set to true, if you want to filter the data based on metric (lower filtering for n_counts and n_genes, and higher filtering for percent_mito and percent_ribo)
   - record_path: path for recording filtered cells CSVs (keep it None if not needed) (default: None) 
    
Requirements:
Pegasus, numpy