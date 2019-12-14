#IMPORTS
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggridges)
library(tools)
source("scripts/mc_functions.R")
source("scripts/fc_plots.R")
source("scripts/readers.R")


#PATHS
data.dir <<- "/broad/hptmp/malperov/data/"
output.dir <<- "/broad/hptmp/malperov/output/"
source.dir <<- output.dir
files.to.save.path <<- "files_to_save.txt"

#FILTERING
cells.filter <<- 3
features.filter <<- 100

#OTHER
save.res.1 <<- TRUE #saveRDS of res 1 automatically

set.seed(29)