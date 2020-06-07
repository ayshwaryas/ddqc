project <<- commandArgs(trailingOnly = TRUE)[1]
task.id <<- as.numeric(Sys.getenv("SGE_TASK_ID")) - 1
source("scripts/readers.R")
source("scripts/filters.R")
source("scripts/mc_functions.R")
source("scripts/settings.R")
source("scripts/local_settings.R")

parse.task.id()
if (method == "none") {
  tiss <- AutoReader(project, cells.filter, features.filter, tasks.per.tiss)
  info.msg <- paste0("task.id:", task.id, " - tissue:", tissue, " res:", res, " mehtod:", method, " param:", param, 
                     " project:", project, " do.counts:", do.counts, " do.genes:", do.genes, " do.mito:", do.mito,
                     " do.ribo:", do.ribo)
  message(paste0("Starting ", info.msg))
  
  task.directory <- paste0(res, "-", method, "-", param)
  task.name <<- paste0(tissue, "-", task.directory)
  robjs.dir <<- paste0(output.dir, "robjs/", project, "/") #directory for saving R objects
  results.dir <<- paste0(output.dir, project, "/", tissue, "/", task.directory, "/") #directory for saving all other output
  
  message("Recording Additional Stats")
  #record all cells with QC stats into csv for summary plots
  t <- tibble("tissue" = tissue, "cell" = colnames(tiss), "nCount_RNA" = tiss$nCount_RNA, "nFeature_RNA" = tiss$nFeature_RNA, "percent.mt" = tiss$percent.mt, "percent.rb" = tiss$percent.rb)
  write.table(t, paste0(results.dir, "../../stats_summary.csv"), sep=",", append=TRUE, col.names=FALSE)
}
