project <<- commandArgs(trailingOnly = TRUE)[1]
task.id <<- as.numeric(Sys.getenv("SGE_TASK_ID")) - 1
source("scripts/readers.R")
source("scripts/mc_functions.R")
source("scripts/fc_plots.R")
source("scripts/settings.R")
source("scripts/local_settings.R")

tiss <- AutoReader(project, cells.filter, features.filter, tasks.per.tiss)
FCPlotsMain()