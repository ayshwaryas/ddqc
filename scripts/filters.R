#filter methods


filterCounts <- function(obj, res, method, threshold) {
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    co <- 0
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$nCount_RNA)[[2]]
      q3 <- summary(cluster$nCount_RNA)[[5]]
      co <- q1 - 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      co <- median(cluster$nCount_RNA) - threshold * mad(cluster$nCount_RNA)
    }
    if (method == "z_score") {
      co <-  mean(cluster$nCount_RNA) - threshold * sd(cluster$nCount_RNA)
    }
    if (method == "cutoff") {
      co <- threshold
    }
    tmp_qc.pass <- (cluster$nCount_RNA >= co)
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  obj[["qc.pass"]] <- qc.pass
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "nCount_RNA" = obj$nCount_RNA[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_counts.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}


filterGenes <- function(obj, res, method, threshold) {
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    lower.co <- 0
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$nFeature_RNA)[[2]]
      q3 <- summary(cluster$nFeature_RNA)[[5]]
      lower.co <- q1 - 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      lower.co <- median(cluster$nFeature_RNA) - threshold * mad(cluster$nFeature_RNA)
    }
    if (method == "z_score") {
      lower.co <- mean(cluster$nFeature_RNA) - threshold * sd(cluster$nFeature_RNA)
    }
    if (method == "cutoff") {
      lower.co <- threshold
    }
    tmp_qc.pass <- (cluster$nFeature_RNA >= lower.co) 
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  obj[["qc.pass"]] <- qc.pass
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "nFeature_RNA" = obj$nFeature_RNA[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_genes.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}


filterMito <- function(obj, res, method, threshold) {
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    co <- 0
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$percent.mt)[[2]]
      q3 <- summary(cluster$percent.mt)[[5]]
      co <- q3 + 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      co <- median(cluster$percent.mt) + threshold * mad(cluster$percent.mt)
    }
    if (method == "z_score") {
      co <- mean(cluster$percent.mt) + threshold * sd(cluster$percent.mt)
    }
    if (method == "cutoff") {
      co <- threshold
    }
    tmp_qc.pass <- (cluster$percent.mt <= co)
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  obj[["qc.pass"]] <- qc.pass
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "percent.mt" = obj$percent.mt[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_mito.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}


filterRibo <- function(obj, res, method, threshold) {
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    co <- 0
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$percent.rb)[[2]]
      q3 <- summary(cluster$percent.rb)[[5]]
      co <- q3 + 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      co <- median(cluster$percent.rb) + threshold * mad(cluster$percent.rb)
    }
    if (method == "z_score") {
      co <-  mean(cluster$percent.rb) + threshold * sd(cluster$percent.rb)
    }
    if (method == "cutoff") {
      co <- threshold
    }
    tmp_qc.pass <- (cluster$percent.rb <= co)
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  obj[["qc.pass"]] <- qc.pass
  tryCatch({ #record filtered cells
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "percent.rb" = obj$percent.rb[filtered.cells])
    write.csv(d, paste0(results.dir, "!filtered_ribo.csv"))
  }, error = function(e) NULL) #in case none of the cells were filtered
  return(qc.pass)
}


filterCells <- function(obj, method, threshold, do.counts, do.genes, do.mito, do.ribo) {
  if (method == "none") {
    return(obj)
  }
  message("Filtering Cells")
  
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
  if (do.counts) {
    if (method == "cutoff") {
      obj$counts.qc.pass <- filterCounts(obj, res, method, 0)
    }
    else {
      obj$counts.qc.pass <- filterCounts(obj, res, method, threshold)
    } 
  }
  else {
    obj$counts.qc.pass <- TRUE
  }
  
  if (do.genes) {
    if (method == "cutoff") {
      obj$genes.qc.pass <- filterGenes(obj, res, method, 200)
    }
    else {
      obj$genes.qc.pass <- filterGenes(obj, res, method, threshold)
    } 
  }
  else {
    obj$genes.qc.pass <- TRUE
  }
  
  if (do.mito) {
    if (method == "cutoff") {
      obj$mito.qc.pass <- filterMito(obj, res, method, threshold)
    }
    else {
      obj$mito.qc.pass <- filterMito(obj, res, method, threshold)
    } 
  }
  else {
    obj$mito.qc.pass <- TRUE
  }
  
  if (do.ribo) {
    if (method == "cutoff") {
      obj$ribo.qc.pass <- filterRibo(obj, res, method, 100)
    }
    else {
      obj$ribo.qc.pass <- filterRibo(obj, res, method, threshold)
    } 
  }
  else {
    obj$ribo.qc.pass <- TRUE
  }
  
  obj <- subset(obj, counts.qc.pass == TRUE & genes.qc.pass == TRUE & mito.qc.pass == TRUE & ribo.qc.pass == TRUE)
  return(obj)
}