source("plotting.R")

PATH <- "/Volumes/easystore/primes_storage/"
DATA.DIR <- paste0(PATH, "output_pg/")

tissues <- c("human_other/Heart_Circulation",
             "human_other/Adipose",
             "human_other/Kidney2",
             "human_other/Liver",
             "human_other/Krasnow_Lung",
             "PanglaoDB/Bone_Marrow",
             "PanglaoDB/Mammary_Gland",
             "PanglaoDB/Pancreatic_Islets",
             "PanglaoDB/Substantia_Nigra",
             "PanglaoDB/Testis",
             "tabula_muris_smartseq2/Bone_Marrow",
             "tabula_muris_smartseq2/Cerebellum",
             "tabula_muris_smartseq2/Colon",
             "tabula_muris_smartseq2/Heart_and_Aorta",
             "tabula_muris/Bladder",
             "tabula_muris/Heart_and_Aorta",
             "tabula_muris/Lung",
             "tabula_muris/Mammary_Gland",
             "tabula_muris/Tongue",
             "tabula_muris/Trachea")

for (tissue in tissues) {
  task.name <<- gsub("/", "-", tissue)
  source.dir <<- paste0(DATA.DIR, tissue, "/1.4-none-0/")
  results.dir <<- paste0(PATH, "figure2_plots/", task.name, "/")
  
  dir.create(paste0(PATH, "figure2_plots/"), showWarnings = FALSE)
  dir.create(paste0(results.dir), showWarnings = FALSE)
  
  file.copy(paste0(source.dir, "!cells.csv"), paste0(results.dir, "!cells.csv"))
  file.copy(paste0(source.dir, "!clusters.csv"), paste0(results.dir, "!clusters.csv"))
  file.copy(paste0(source.dir, "!markers.csv"), paste0(results.dir, "!markers.csv"))
  
  #cells <- read.csv(paste0(source.dir, "!cells.csv"))
  #cells$louvain_labels <- as.factor(cells$louvain_labels)
  #clusters <- read.csv(paste0(source.dir, "!clusters.csv"))
  
  #generatePlots(cells, clusters$cell_type, no.ct.labels=TRUE)
}