project <<- commandArgs(trailingOnly = TRUE)[1]
task.id <<- as.numeric(Sys.getenv("SGE_TASK_ID")) - 1
source("scripts/readers.R")
source("scripts/mc_functions.R")
source("scripts/settings.R")
source("scripts/local_settings.R")

tasks.per.res <- tasks.per.tiss %/% 2 #how many different methods per one resolution
res <<- 0.5 * (1 + ((task.id %% tasks.per.tiss) %/% tasks.per.res)) #clustering resolution
method <<- switch(task.id %% tasks.per.res + 1, "none", "cutoff", "z_score", "cutoff", "cutoff", "outlier", "percentile", "mad") #filtering method
param <<- switch(task.id %% tasks.per.res + 1, 0, 80, 2, 5, 10, 0, 95, 2) #filtering parameter

if (res == 0.5 && method == "none") { 
  message("Recording Additional Stats")
  filename <- paste0(res, "-", method, "-", param)
  robjs.dir <<- paste0(output.dir, "robjs/", project, "/") #directory for saving R objects
  results.dir <<- paste0(output.dir, project, "/", tissue, "/", filename, "/") #directory for saving all other output
  tiss <- AutoReader(project, cells.filter, features.filter, tasks.per.tiss)
  #record all cells with QC stats into csv for summary plots
  t <- tibble("tissue" = tissue, "nCount_RNA" = tiss$nCount_RNA, "nFeature_RNA" = tiss$nFeature_RNA, "percent.mt" = tiss$percent.mt, "percent.rb" = tiss$percent.rb)
  write.table(t, paste0(results.dir, "../../stats_summary.csv"), sep=",", append=TRUE, col.names=FALSE)
}
