source("plots/tmp_utils.R")
source("plots/plotting.R")


task.name <<- commandArgs(trailingOnly = TRUE)[1]
results.dir <<- commandArgs(trailingOnly = TRUE)[2]
message("Starting R script to generate results")

results.dir <<- paste0(output.dir, project, "/", tissue, "/filtered_cells_plots/", metric, "/") #directory for saving all other output

tiss <- read.csv(paste0(results.dir, "!cells.csv"))

generatePlots(tiss, task.name, clusters$cell.type)
generateFCPlots(tiss, clusters)

