#IMPORTS
library(dplyr)
if (!exists("data.from.pg")){
  library(Seurat)
}
library(ggplot2)
library(cowplot)
library(ggridges)
library(tools)


#PATHS
data.dir <<- "/ahg/regevdata/projects/scqc/data/"
if (! exists("data.from.pg")) {
  output.dir <<- "/ahg/regevdata/projects/scqc/output/"
} else {
  output.dir <<- "/ahg/regevdata/projects/scqc/output_pg/"
}
source.dir.prefix <<- output.dir

#FILTERING
cells.filter <<- 3
features.filter <<- 100
do.counts <<- TRUE
do.genes <<- TRUE
do.mito <<- TRUE
do.ribo <<- FALSE

#OTHER
save.res.1 <<- FALSE #saveRDS of res 1 automatically

set.seed(29)