source("scripts/readers.R")
source("scripts/mc_functions.R")
source("scripts/fc_plots.R")
source("scripts/settings.R")
source("scripts/local_settings.R")
tasks.per.tiss <<- 1


project <- "mc_tm"
tissue <- "Lung"


readFilterCsvMethod <- function(method, all.cells) {
  filtered.cells.all <- NULL
  line <- method
  for (metric in c("counts", "genes", "mito", "ribo")) {
    filtered.cells <- tryCatch({
      as.character(read.csv(paste0(source.dir, res, "-", method, "/!filtered_", metric, ".csv"))[["barcodekey"]])}, 
      error = function(e) {warning(paste(method, metric, "filtered cells not found"))})
    filtered.cells.all <- union(filtered.cells.all, filtered.cells)
    line <- paste(line, length(all.cells) - length(filtered.cells), paste0(
      round((length(all.cells) - length(filtered.cells)) / length(all.cells) * 100, 1), "%"), sep=",")
  }
  line <- paste(line, length(all.cells) - length(filtered.cells.all), paste0(
    round((length(all.cells) - length(filtered.cells.all)) / length(all.cells) * 100, 1), "%"), sep=",")
  write(line, file=paste0(source.dir, "fc_table.csv"), append=TRUE)
  return(filtered.cells)
}



res <<- 1.4 #0.5 * (1 + (task.id %% tasks.per.tiss) %/% tasks.per.res) #clustering resolution
source.dir <<- paste0(source.dir.prefix, project, "/", tissue, "/") #directory where csv with filtered cells are located

all.cells = as.character(read.csv(paste0(source.dir, res, "-none-0/!cells.csv"))$barcodekey)

write("method,counts cells,counts %,genes cells,genes %,mito cells,mito %,ribo cells,ribo %,all cells,all %", file=paste0(source.dir, "fc_table.csv"))

#cutoff5 <- setdiff(all.cells, readFilterCsvMethod("cutoff-5", all.cells))
cutoff10 <- setdiff(all.cells, readFilterCsvMethod("cutoff-10", all.cells))
#zscore2 <- setdiff(all.cells, readFilterCsvMethod("z_score-2", all.cells))
mad <- setdiff(all.cells, readFilterCsvMethod("mad-2", all.cells))
outlier <- setdiff(all.cells, readFilterCsvMethod("outlier-0", all.cells))
   

