#filter methods


#//////// NCounts ///////////
adaptiveFilter.counts <- function(obj, res, method, param) {
  original <- obj #make a copy of the object 
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    co <- 0
    if (method == "percentile") { #percentile filtering
      co <-  quantile(cluster$nCount_RNA, probs = (1 - param / 100))
    }
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$nCount_RNA)[[2]]
      q3 <- summary(cluster$nCount_RNA)[[5]]
      co <- q1 - 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      co <- median(obj$nCount_RNA) - param * mad(cluster$nCount_RNA)
    }
    if (method == "z_score") {
      co <-  mean(cluster$nCount_RNA) - param * sd(cluster$nCount_RNA)
    }
    tmp_qc.pass <- (cluster$nCount_RNA >= co)
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  original[["qc.pass"]] <- qc.pass
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(original, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = original$annotations[filtered.cells], "nCount_RNA" = original$nCount_RNA[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_counts.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}

cutoffFilter.counts <- function(obj, param) { #filter cells with %counts > param
  obj[["qc.pass"]] <- (obj$nCount_RNA >= param)
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "nCount_RNA" = obj$nCount_RNA[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_counts.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}


#//////// NGene ///////////
adaptiveFilter.genes <- function(obj, res, method, param) {
  original <- obj #make a copy of the object 
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    lower.co <- 0
    upper.co <- 0
    if (method == "percentile") { #percentile filtering
      lower.co <- quantile(cluster$nFeature_RNA, probs = (100 - param) / 100)
      upper.co <- quantile(cluster$nFeature_RNA, probs = param / 100)
    }
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$nFeature_RNA)[[2]]
      q3 <- summary(cluster$nFeature_RNA)[[5]]
      lower.co <- q1 - 1.5 * (q3 - q1)
      upper.co <- q3 + 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      lower.co <- median(cluster$nFeature_RNA) - param * mad(cluster$nFeature_RNA)
      upper.co <- median(cluster$nFeature_RNA) + param * mad(cluster$nFeature_RNA)
    }
    if (method == "z_score") {
      lower.co <- mean(cluster$nFeature_RNA) - param * sd(cluster$nFeature_RNA)
      upper.co <- mean(cluster$nFeature_RNA) + param * sd(cluster$nFeature_RNA)
    }
    tmp_qc.pass <- (cluster$nFeature_RNA >= lower.co & cluster$nFeature_RNA <= upper.co) 
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  original[["qc.pass"]] <- qc.pass
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(original, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = original$annotations[filtered.cells], "nFeature_RNA" = original$nFeature_RNA[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_genes.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}

cutoffFilter.genes <- function(obj, param1, param2) { #filter cells by nFeature_RNA
  obj[["qc.pass"]] <- (obj$nFeature_RNA >= param1 & obj$nFeature_RNA <= param2)
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "nFeature_RNA" = obj$nFeature_RNA[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_genes.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}

#//////// %Mito ///////////
adaptiveFilter.mito <- function(obj, res, method, param) {
  original <- obj #make a copy of the object 
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    co <- 0
    if (method == "percentile") { #percentile filtering
      co <-  quantile(cluster$percent.mt, probs = param / 100)
    }
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$percent.mt)[[2]]
      q3 <- summary(cluster$percent.mt)[[5]]
      co <- q3 + 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      co <- median(obj$percent.mt) + param * mad(cluster$percent.mt)
    }
    if (method == "z_score") {
      co <-  mean(cluster$percent.mt) + param * sd(cluster$percent.mt)
    }
    tmp_qc.pass <- (cluster$percent.mt <= co)
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  original[["qc.pass"]] <- qc.pass
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(original, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = original$annotations[filtered.cells], "percent.mt" = original$percent.mt[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_mito.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}

cutoffFilter.mito <- function(obj, param) { #filter cells with %mito > param
  obj[["qc.pass"]] <- (obj$percent.mt <= param)
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "percent.mt" = obj$percent.mt[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_mito.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}


#//////// %Ribo ///////////
adaptiveFilter.ribo <- function(obj, res, method, param) {
  original <- obj #make a copy of the object 
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    co <- 0
    if (method == "percentile") { #percentile filtering
      co <-  quantile(cluster$percent.rb, probs = param / 100)
    }
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$percent.rb)[[2]]
      q3 <- summary(cluster$percent.rb)[[5]]
      co <- q3 + 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      co <- median(obj$percent.rb) + param * mad(cluster$percent.rb)
    }
    if (method == "z_score") {
      co <-  mean(cluster$percent.rb) + param * sd(cluster$percent.rb)
    }
    tmp_qc.pass <- (cluster$percent.rb <= co)
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  original[["qc.pass"]] <- qc.pass
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(original, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = original$annotations[filtered.cells], "percent.rb" = original$percent.rb[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_ribo.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}

cutoffFilter.ribo <- function(obj, param) { #filter cells with %ribo > param
  obj[["qc.pass"]] <- (obj$percent.rb <= param)
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "percent.rb" = obj$percent.rb[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_ribo.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}


filterCells <- function(obj, method, param, do.counts, do.genes, do.mito, do.ribo) {
  message("Filtering Cells")
  
  obj[["qc.pass"]] <- TRUE
  
  if (do.counts) {
    if (method == "outlier" || method == "percentile" || method == "mad" || method == "z_score") {
      obj$qc.pass <- obj$qc.pass & adaptiveFilter.counts(obj, res, method, param)
    } 
    if (method == "cutoff" && param < 20) {
      obj$qc.pass <- obj$qc.pass & cutoffFilter.counts(obj, 500)
    }
  }
  
  if (do.genes) {
    if (method == "outlier" || method == "percentile" || method == "mad" || method == "z_score") {
      obj$qc.pass <- obj$qc.pass & adaptiveFilter.genes(obj, res, method, param)
    } 
    if (method == "cutoff" && param < 20) {
      obj$qc.pass <- obj$qc.pass & cutoffFilter.genes(obj, 200, 2500)
    }
  }
  
  if (do.mito) {
    if (method == "outlier" || method == "percentile" || method == "mad" || method == "z_score") {
      obj$qc.pass <- obj$qc.pass & adaptiveFilter.mito(obj, res, method, param)
    } 
    if (method == "cutoff") {
      obj$qc.pass <- obj$qc.pass & cutoffFilter.mito(obj, param)
    }
  }
  
  if (do.ribo) {
    if (method == "outlier" || method == "percentile" || method == "mad" || method == "z_score") {
      obj$qc.pass <- obj$qc.pass & adaptiveFilter.ribo(obj, res, method, param)
    } 
    if (method == "cutoff") {
      obj$qc.pass <- obj$qc.pass & cutoffFilter.ribo(obj, 100)
    }
  }
  
  obj <- subset(obj, qc.pass == TRUE)
  return(obj)
}