#IMPORTS
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggridges)
library(tools)


#PATHS
data.dir <<- "/broad/hptmp/malperov/data/"
output.dir <<- "/ahg/regevdata/projects/scqc/output/"
source.dir <<- output.dir

#FILTERING
cells.filter <<- 3
features.filter <<- 100

#OTHER
save.res.1 <<- FALSE #saveRDS of res 1 automatically

set.seed(29)