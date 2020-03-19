#filter methods


trueMAD <- function(x) {
  return(median(abs(x - median(x))))
}


filterCounts <- function(obj, res, method, param) {
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
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
      co <- median(cluster$nCount_RNA) - param * trueMAD(cluster$nCount_RNA)
    }
    if (method == "z_score") {
      co <-  mean(cluster$nCount_RNA) - param * sd(cluster$nCount_RNA)
    }
    if (method == "cutoff") {
      co <- param
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


filterGenes <- function(obj, res, method, param) {
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) { #for each cluster calculate which cells passed QC
    cluster <- subset(obj, idents = cl)
    lower.co <- 0
    upper.co <- 0
    if (method == "outlier") { #outlier filtering (> Q3 + 1.5IQR)
      q1 <- summary(cluster$nFeature_RNA)[[2]]
      q3 <- summary(cluster$nFeature_RNA)[[5]]
      lower.co <- q1 - 1.5 * (q3 - q1)
      upper.co <- q3 + 1.5 * (q3 - q1)
    }
    if (method == "mad") { #Median absolute deviations filtering 
      lower.co <- median(cluster$nFeature_RNA) - param * trueMAD(cluster$nFeature_RNA)
      upper.co <- median(cluster$nFeature_RNA) + param * trueMAD(cluster$nFeature_RNA)
    }
    if (method == "z_score") {
      lower.co <- mean(cluster$nFeature_RNA) - param * sd(cluster$nFeature_RNA)
      upper.co <- mean(cluster$nFeature_RNA) + param * sd(cluster$nFeature_RNA)
    }
    if (method == "cutoff") {
      lower.co <- param
      upper.co <- param * 10
    }
    tmp_qc.pass <- (cluster$nFeature_RNA >= lower.co & cluster$nFeature_RNA <= upper.co) 
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


filterMito <- function(obj, res, method, param) {
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
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
      co <- median(cluster$percent.mt) + param * trueMAD(cluster$percent.mt)
    }
    if (method == "z_score") {
      co <-  mean(cluster$percent.mt) + param * sd(cluster$percent.mt)
    }
    if (method == "cutoff") {
      co <- param
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


filterRibo <- function(obj, res, method, param) {
  obj <- clusterize(obj, res, compute.reductions = FALSE, compute.markers = FALSE) #cluster the data
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
      co <- median(cluster$percent.rb) + param * trueMAD(cluster$percent.rb)
    }
    if (method == "z_score") {
      co <-  mean(cluster$percent.rb) + param * sd(cluster$percent.rb)
    }
    if (method == "cutoff") {
      co <- param
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


filterCells <- function(obj, method, param, do.counts, do.genes, do.mito, do.ribo) {
  message("Filtering Cells")
  
  obj[["qc.pass"]] <- TRUE
  
  if (do.counts) {
    if (method == "cutoff") {
      obj$qc.pass <- obj$qc.pass & filterCounts(obj, res, method, 500)
    }
    else {
      obj$qc.pass <- obj$qc.pass & filterCounts(obj, res, method, param)
    } 
  }
  
  if (do.genes) {
    if (method == "cutoff") {
      obj$qc.pass <- obj$qc.pass & filterGenes(obj, res, method, 250)
    }
    else {
      obj$qc.pass <- obj$qc.pass & filterGenes(obj, res, method, param)
    } 
  }
  
  if (do.mito) {
    if (method == "cutoff") {
      obj$qc.pass <- obj$qc.pass & filterMito(obj, res, method, param)
    }
    else {
      obj$qc.pass <- obj$qc.pass & filterMito(obj, res, method, param)
    } 
  }
  
  if (do.ribo) {
    if (method == "cutoff") {
      obj$qc.pass <- obj$qc.pass & filterRibo(obj, res, method, 100)
    }
    else {
      obj$qc.pass <- obj$qc.pass & filterRibo(obj, res, method, param)
    } 
  }
  
  obj <- subset(obj, qc.pass == TRUE)
  return(obj)
}