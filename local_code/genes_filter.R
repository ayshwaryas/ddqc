library(dplyr)
library(Seurat)

outlierFilter <- function(obj, res, param, filename) {
  obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x = obj)
  obj <- ScaleData(object = obj, features = all.genes)
  obj <- RunPCA(object = obj, features = VariableFeatures(object = obj))
  obj <- FindNeighbors(object = obj, dims = 1:10)
  obj <- FindClusters(object = obj, resolution = res)
  qc.pass <- NULL
  for (cl in levels(obj$seurat_clusters)) {
    cluster <- subset(obj, idents = cl)
    co
    if (param == 95) {
      co <-  quantile(cluster$nFeature_RNA, probs = 0.05)
    }
    else {
      q1 <- summary(cluster$nFeature_RNA)[[2]]
      q3 <- summary(cluster$nFeature_RNA)[[5]]
      co <- q1 - 1.5 * (q3 - q1)
    }
    tmp_qc.pass <- (cluster$nFeature_RNA >= co)
    names(tmp_qc.pass) <- colnames(cluster$RNA)
    qc.pass <- c(qc.pass, tmp_qc.pass)
  }
  obj[["qc.pass"]] <- qc.pass
  tryCatch({
    filtered.cells <- colnames(subset(obj, qc.pass == FALSE))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "nFeature_RNA" = obj$nFeature_RNA[filtered.cells])
    write.csv(d, paste0(output.dir, filename, "/!filtered.csv"))
  }, error = function(e) NULL)
  obj <- subset(obj, qc.pass == TRUE)
  return(obj)
}

zScoreFilter <- function(obj, res, param, filename) {
  obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
  # linear dimension reduction
  obj <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x = obj)
  obj <- ScaleData(object = obj, features = all.genes)
  obj <- RunPCA(object = obj, features = VariableFeatures(object = obj))
  obj <- FindNeighbors(object = obj, dims = 1:10)
  obj <- FindClusters(object = obj, resolution = res)
  z_scores <- NULL
  for (cl in levels(obj$seurat_clusters)) {
    cluster <- subset(obj, idents = cl)
    tmp_z_scores <- (cluster$nFeature_RNA - mean(cluster$nFeature_RNA)) / sd(cluster$nFeature_RNA)
    names(tmp_z_scores) <- colnames(cluster$RNA)
    z_scores <- c(z_scores, tmp_z_scores)
  }
  obj[["z_scores"]] <- z_scores
  tryCatch({
    filtered.cells <- colnames(subset(obj, z_scores > param))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "nFeature_RNA" = obj$nFeature_RNA[filtered.cells])
    write.csv(d, paste0(output.dir, filename, "/!filtered.csv"))
  }, error = function(e) NULL)
  obj <- subset(obj, z_scores >= -param)
  return(obj)
}

cutoffFilter <- function(obj, param, filename) {
  tryCatch({
    filtered.cells <- colnames(subset(obj, nFeature_RNA > param))
    d <- tibble("cell" = filtered.cells, "annotation" = obj$annotations[filtered.cells], "nFeature_RNA" = obj$nFeature_RNA[filtered.cells])
    write.csv(d, paste0(output.dir, filename, "/!filtered.csv"))
  }, error = function(e) NULL)
  obj <- subset(obj, nFeature_RNA >= param)
  return(obj)
}

