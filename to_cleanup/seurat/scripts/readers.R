ReadTabulaMuris <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "Bladder", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland", "Marrow", "Spleen", "Thymus", "Tongue", "Trachea")
  is.human <<- FALSE
  data.path <<- paste0(data.dir, "tabula_muris/droplet/")
  files <- list.files(path=data.path, pattern=tissue)
  tiss <- NULL
  if (length(files) > 1) {
    objs <- NULL
    filenames <- NULL
    for (file in files) {
      objs <- c(objs, CreateSeuratObject(counts = Read10X(data.dir=paste0(data.path, file)), project = (strsplit(file, "-")[[1]][[2]]), min.cells = cells.filter, min.features = features.filter))
      filenames <- c(filenames, (strsplit(file, "-")[[1]][[2]]))
    }
    tiss <- merge(x = objs[[1]], y = objs[2:length(objs)], add.cell.ids = filenames, project = tissue)
  } else {
    file <- files[1]
    tiss <- CreateSeuratObject(counts = Read10X(data.dir = paste0(data.path, file)), project = (strsplit(file, "-")[[1]][[2]]), min.cells = cells.filter, min.features = features.filter)
    tiss <- RenameCells(tiss, add.cell.id = (strsplit(file, "-")[[1]][[2]]))
  }
  
  tiss <- RenameCells(tiss, new.names = sapply(strsplit(colnames(tiss), "-"), first))
  annotations <- read.csv(file= paste0(data.path, "annotations_droplet.csv"), header=TRUE, sep=",")
  annotations.cell.type <- as.character(annotations[["cell_ontology_class"]])
  names(annotations.cell.type) <- as.character(annotations[["cell"]])
  tiss[["annotations"]] <- annotations.cell.type
  tiss$annotations[is.na(tiss$annotations)] <- "Unknown"
  return(tiss)
}

ReadTabulaSenis24 <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "Bladder", "Bm", "Gat",  "Heart", "Hepatocyctes", "Kidney", "Mat", "Muscle", "Pancrease", "Scat", "Spleen", "Thymus", "Tongue")
  is.human <<- FALSE
  data.path <<- paste0(data.dir, "tabula_senis/10X/24_month/")
  files <- list.files(path=data.path, pattern=paste0("M_", tissue), ignore.case = TRUE)
  tiss <- NULL
  if (length(files) > 1) {
    objs <- NULL
    filenames <- NULL
    for (file in files) {
      objs <- c(objs, CreateSeuratObject(counts = Read10X_h5(file = paste0(data.path, file, "/raw_gene_bc_matrices_h5.h5")), project = file, min.cells = cells.filter, min.features = features.filter))
      filenames <- c(filenames, file)
    }
    tiss <- merge(x = objs[[1]], y = objs[2:length(objs)], add.cell.ids = filenames, project = tissue)
  } else {
    file <- files[1]
    tiss <- CreateSeuratObject(counts = Read10X_h5(file = paste0(data.path, file, "/raw_gene_bc_matrices_h5.h5")), project = file, min.cells = cells.filter, min.features = features.filter)
    tiss <- RenameCells(tiss, add.cell.id = file)
  }
  
  annotations <- read.csv(file= paste0(data.path, "annotations.csv"), header=TRUE, sep=",")
  annotations.cell.type <- as.character(annotations[["cell_ontology_class_reannotated"]])
  names(annotations.cell.type) <- as.character(paste0(annotations[["cell"]], "-1"))
  tiss[["annotations"]] <- annotations.cell.type
  tiss$annotations[is.na(tiss$annotations)] <- "Unknown"
  tiss$annotations[tiss$annotations == "nan"] <- "Unknown"
  return(tiss)
}

ReadTabulaSenis30 <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "Fat", "Heart_and_Aorta", "Kidney", "Large_Intestine", "Limb_Muscle", "Liver", "Lung", "Marrow", "Pancreas", "Spleen")
  is.human <<- FALSE
  data.path <<- paste0(data.dir, "tabula_senis/10X/30_month/")
  files <- list.files(path=data.path, pattern=tissue)
  tiss <- NULL
  if (length(files) > 1) {
    objs <- NULL
    filenames <- NULL
    for (file in files) {
      objs <- c(objs, CreateSeuratObject(counts = Read10X_h5(file = paste0(data.path, file, "/raw_gene_bc_matrices_h5.h5")), project = (strsplit(file, "-")[[1]][[2]]), min.cells = cells.filter, min.features = features.filter))
      filenames <- c(filenames, (strsplit(file, "-")[[1]][[2]]))
    }
    tiss <- merge(x = objs[[1]], y = objs[2:length(objs)], add.cell.ids = filenames, project = tissue)
  } else {
    file <- files[1]
    tiss <- CreateSeuratObject(counts = Read10X_h5(file = paste0(data.path, file, "/raw_gene_bc_matrices_h5.h5")), project = (strsplit(file, "-")[[1]][[2]]), min.cells = cells.filter, min.features = features.filter)
    tiss <- RenameCells(tiss, add.cell.id = (strsplit(file, "-")[[1]][[2]]))
  }
  
  annotations <- read.csv(file= paste0(data.path, "annotations.csv"), header=TRUE, sep=",")
  annotations.cell.type <- as.character(annotations[["cell_ontology_class_reannotated"]])
  names(annotations.cell.type) <- as.character(paste0(annotations[["cell"]], "-1"))
  tiss[["annotations"]] <- annotations.cell.type
  tiss$annotations[is.na(tiss$annotations)] <- "Unknown"
  tiss$annotations[tiss$annotations == "nan"] <- "Unknown"
  return(tiss)
}

ReadMouseCellAtlas <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "Bladder", "BoneMarrow", "Brain", "Kidney", "Liver", "Lung", "MammaryGland.Involution", "MammaryGland.Lactation", "MammaryGland.Pregnancy", "MammaryGland.Virgin", "MesenchymalStemCells", "Muscle", "Ovary", "Pancreas", "PeripheralBlood", "Placenta", "Prostate", "SmallIntestine", "Spleen", "Stomach", "Testis", "Thymus", "TrophoblastStemCells", "Uterus")
  is.human <<- FALSE
  data.path <<- paste0(data.dir, "mca/")
  tiss <- NULL
  files <- list.files(path=data.path, pattern=paste0("^", tissue))
  if (length(files) > 1) {
    objs <- NULL
    for (file in files) {
      objs <- c(objs, CreateSeuratObject(counts = read.table(paste0(data.path, file)), min.cells = cells.filter, min.features = features.filter))
    }
    tiss <- merge(x = objs[[1]], y = objs[2:length(objs)])
    rm(objs)
  } else {
    tiss <- CreateSeuratObject(counts = read.table(paste0(data.path, files[[1]])), min.cells = cells.filter, min.features = features.filter)
  }
  
  annotations <- read.csv(file= paste0(data.path, "MCA_CellAssignments.csv"), header=TRUE, sep=",")
  annotations.annotation <- as.character(annotations[["Annotation"]])
  names(annotations.annotation) <- as.character(annotations[["Cell.name"]])
  tiss[["annotations"]] <- sapply(strsplit(annotations.annotation, "(\\()|(_)"), first)
  tiss$annotations[is.na(tiss$annotations)] <- "Unknown"
  return(tiss)
}

ReadEBI <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "E-ENAD-21-breast_epithelial1", "E-ENAD-27-islet", "E-GEOD-81547-pancreas", "E-GEOD-81608-islet2", "E-GEOD-83139-pancreatic_endocrine", "E-GEOD-86618-lung_epithelial", "E-GEOD-89232-dendritic_cells", "E-GEOD-130148-lung", "E-MTAB-6386-B_cells", "E-MTAB-6653-lung_carcinomas", "E-MTAB-6701-fetal-maternal_interface", "E-CURD-6-bone_marrow", "E-MTAB-7316-retina")
  is.human <<- TRUE
  data.path <<- paste0(data.dir, "ebi/")
  file <- list.files(path=data.path, pattern=tissue)[1]
  tiss <- CreateSeuratObject(counts = Read10X(data.dir = paste0(data.path, file)), project = tissue, min.cells = cells.filter, min.features = features.filter)
  annotations <- read.csv(file=paste0(data.path, file, "/annotations.tsv"), sep="\t", header=TRUE)
  ind <- which(grepl("cell.*type", colnames(annotations)), TRUE)
  if (length(ind) > 0) {
    ind <- ind[[1]]
    annotations.cell.type <- as.character(annotations[[ind]])
    names(annotations.cell.type) <- as.character(annotations[["Assay"]])
    tiss[["annotations"]] <- annotations.cell.type
    tiss$annotations[is.na(tiss$annotations)] <- "not available"
  } else {
    tiss[["annotations"]] <- "Unknown"
  }
  return(tiss)
}

ReadEBI_TM <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "Adipose", "Bladder", "Bone_Marrow", "Cerebellum", "Cerebral_Cortex", "Colon", "Diaphragm", "Fat", "Heart_and_Aorta", "Hippocampus", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland", "Pancreas", "Skin", "Spleen", "Striatum", "Thymus", "Tongue", "Trachea")
  is.human <<- FALSE
  data.path <<- paste0(data.dir, "ebi_tm/")
  file <- paste0(data.path, tissue, ".rds")
  tiss.unfiltered <- readRDS(file)
  annotations <- tiss.unfiltered$annotations
  tiss <- CreateSeuratObject(tiss.unfiltered$RNA@counts, min.cells=cells.filter, min.features=features.filter)
  tiss[["annotations"]] <- annotations
  return(tiss)
}

ReadOther <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "bipolar_mouse", "colon-epi_human", "colon-fib_human", "colon-imm_human", "kidney_human", "retina_human", "retina_mouse")
  is.human <<- grepl("human", tissue)
  data.path <<- paste0(data.dir, "other/")
  file <- paste0(data.path, tissue, ".rds")
  tiss.unfiltered <- readRDS(file)
  annotations <- tiss.unfiltered$annotations
  tiss <- CreateSeuratObject(tiss.unfiltered$RNA@counts, min.cells=cells.filter, min.features=features.filter)
  tiss[["annotations"]] <- annotations
  return(tiss)
}

ReadOther10X <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "adipose", "ASD_snRNAseq", "liver", "skin")
  is.human <<- TRUE
  data.path <<- paste0(data.dir, "other/", tissue, "/")
  files <- list.files(path=data.path, pattern = tissue)
  tiss <- NULL
  if (length(files) > 1) {
    objs <- NULL
    filenames <- NULL
    for (file in files) {
      objs <- c(objs, CreateSeuratObject(counts = Read10X(data.dir = paste0(data.path, file)), min.cells = cells.filter, min.features = features.filter))
      filenames <- c(filenames, file)
    }
    tiss <- merge(x = objs[[1]], y = objs[2:length(objs)], project = tissue)
  } else {
    file <- files[1]
    tiss <- CreateSeuratObject(counts = Read10X(data.dir = paste0(data.path, file)), min.cells = cells.filter, min.features = features.filter)
  }
  tiss[["annotations"]] <- "Unknown"
  if (tissue == "adipose" || tissue == "liver") {
    ann <- read.csv(paste0(data.path, "annotations.csv"))
    annotations <- as.character(ann$x)
    names(annotations) <- as.character(ann$X)
    tiss$annotations <- annotations
    tiss$annotations[is.na(tiss$annotations)] <- "Unknown"
  }
  return(tiss)
}

ReadKB_TM <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "Kidney", "Trachea")
  is.human <<- FALSE
  data.path <<- paste0(data.dir, "kb_test/results/", tissue, "/")
  files <- list.files(data.path)
  tiss <- NULL
  if (length(files) > 1) {
    objs <- NULL
    filenames <- NULL
    for (file in files) {
      objs <- c(objs, CreateSeuratObject(readRDS(paste0(data.path, file))$RNA@counts, min.cells=cells.filter, min.features=features.filter))
      filenames <- c(filenames, strsplit(file, "\\.")[[1]][1])
    }
    tiss <- merge(x = objs[[1]], y = objs[2:length(objs)], add.cell.ids = filenames, project = tissue)
  } else {
    file <- files[1]
    tiss <- CreateSeuratObject(readRDS(paste0(data.path, file))$RNA@counts, min.cells=cells.filter, min.features=features.filter)
    tiss <- RenameCells(tiss, add.cell.id = strsplit(file, "\\.")[[1]][1])
  }
  
  annotations <- read.csv(file= paste0(data.dir, "tabula_muris/droplet/annotations_droplet.csv"), header=TRUE, sep=",")
  annotations.cell.type <- as.character(annotations[["cell_ontology_class"]])
  names(annotations.cell.type) <- as.character(annotations[["cell"]])
  tiss[["annotations"]] <- annotations.cell.type
  tiss$annotations[is.na(tiss$annotations)] <- "Unknown"
  return(tiss)
}

ReadPanglaoDB <- function(cells.filter, features.filter, tasks.per.tiss) {
  tissue <<- switch(1 + (task.id %/% tasks.per.tiss), "Lung_Epithelial", "Mammary_Gland", "Pancreatic_Islets", "Prostate", "Testis")
  is.human <<- TRUE
  data.path <<- paste0(data.dir, "PanglaoDB/")
  files <- list.files(path=data.path, pattern=tissue)
  tiss <- NULL
  if (length(files) > 1) {
    objs <- NULL
    filenames <- NULL
    for (file in files) {
      load(paste0(data.path, file))
      gene.names <- rownames(sm)
      gene.names <- sapply(strsplit(gene.names, "_ENSG"), first)
      rownames(sm) <- gene.names
      objs <- c(objs, CreateSeuratObject(counts = sm, project = (strsplit(file, "-")[[1]][[2]]), min.cells = cells.filter, min.features = features.filter))
      filenames <- c(filenames, (strsplit(file, "-")[[1]][[2]]))
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
  return(tiss)
}

AutoReader <- function(dataset, cells.filter, features.filter, tasks.per.tiss) {
  message(paste("Reading Data - dataset:", dataset, "min.cells:", cells.filter, "min.features:", features.filter))
  if (dataset == "tabula_muris" || dataset == "mc_tm" || dataset == "tm") {
    obj <- ReadTabulaMuris(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "tabula_senis24" || dataset == "mc_ts24" || dataset == "ts24") {
    obj <- ReadTabulaSenis24(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "tabula_senis30" || dataset == "mc_ts30"|| dataset == "ts30") {
    obj <- ReadTabulaSenis30(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "mca"|| dataset == "mc_mca") {
    obj <- ReadMouseCellAtlas(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "ebi"|| dataset == "mc_ebi") {
    obj <- ReadEBI(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "other" || dataset == "mc_other") {
    obj <- ReadOther(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "other_10X" || dataset == "mc_other_10X") {
    obj <- ReadOther10X(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "ebi_tm" || dataset=="mc_ebi_tm") {
    obj <- ReadEBI_TM(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "kb_tm" || dataset=="mc_kb_tm") {
    obj <- ReadKB_TM(cells.filter, features.filter, tasks.per.tiss)
  }
  if (dataset == "PanglaoDB" || dataset=="mc_PanglaoDB") {
    obj <- ReadPanglaoDB(cells.filter, features.filter, tasks.per.tiss)
  }
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, features=grep("^MT-", rownames(obj$RNA), ignore.case=TRUE))
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, features=grep("^RP[SL][[:digit:]]", rownames(obj$RNA), ignore.case=TRUE))
  
  obj <- subset(obj, percent.mt <= 80)
  
  return(obj)
}