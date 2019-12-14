project <<- Sys.getenv("R_PROJECT")
task.id <<- as.numeric(Sys.getenv("SGE_TASK_ID")) - 1
source("scripts/settings.R")

tiss <- AutoReader(project, cells.filter, features.filter, tasks.per.tiss)
MCMain(tiss)
