#IMPORTS
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggridges)
library(tools)


#PATHS
data.dir <<- "/broad/hptmp/malperov/data/"
output.dir <<- "/broad/hptmp/malperov/output/"
source.dir <<- output.dir
files.to.save.path <<- "files_to_save.txt"

#FILTERING
cells.filter <<- 3
features.filter <<- 100

#OTHER
save.res.1 <<- FALSE #saveRDS of res 1 automatically

set.seed(29)