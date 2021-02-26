library(Seurat)
library(Matrix)

first <- function(x) {
  return(x[1])
}


data.dir <- "~/Documents/primes_storage/data/PanglaoDB/"
features.filter <- 0 
cells.filter <- 0
for (tissue in c("Pancreatic_Islets", "Prostate", "Testis")) {
  is.human <<- TRUE
  data.path <<- paste0(data.dir, tissue, "/")
  files <- list.files(path=data.path, pattern="SRA")
  tiss <- NULL
  if (length(files) > 1) {
    objs <- NULL
    filenames <- NULL
    for (file in files) {
      load(paste0(data.path, file))
      gene.names <- rownames(sm)
      gene.names <- sapply(strsplit(gene.names, "_ENSG"), first)
      rownames(sm) <- gene.names
      objs <- c(objs, CreateSeuratObject(counts = sm, project = (file), min.cells = cells.filter, min.features = features.filter))
      filenames <- c(filenames, file)
    }
    tiss <- merge(x = objs[[1]], y = objs[2:length(objs)], add.cell.ids = filenames, project = tissue)
  } else {
    file <- files[1]
    load(paste0(data.path, file))
    gene.names <- rownames(sm)
    gene.names <- sapply(strsplit(gene.names, "_ENSG"), first)
    rownames(sm) <- gene.names
    tiss <- CreateSeuratObject(counts = sm, project = (strsplit(file, "-")[[1]][[2]]), min.cells = cells.filter, min.features = features.filter)
    tiss <- RenameCells(tiss, add.cell.id = (strsplit(file, "-")[[1]][[2]]))
  }
  
  tiss$annotations <- "Unknown"
  
  
  sparse.gbm <- Matrix(tiss$RNA@counts, sparse = T )
  
  directory <- paste0(data.dir, tissue, "/")
  dir.create(directory)
  
  writeMM(obj = sparse.gbm, file=paste0(directory, "matrix.mtx"))
  write(x = rownames(tiss$RNA@counts), file=paste0(directory, "genes.tsv"))
  write(x = colnames(tiss$RNA@counts), file=paste0(directory, "barcodes.tsv"))
  
  annotations <- as.data.frame(tiss$annotations)
  colnames(annotations) <- c("celltype")
  write.csv(annotations, file=paste0(directory, "annotations.csv"))
  
}


