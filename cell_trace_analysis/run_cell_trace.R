library(optparse)

option_list <- list(
  make_option(c('-s', '--seurat_Obj'), type = 'character', help = 'seurat Object'),
  make_option(c('-p','--species'), type = 'character', 
              help = 'species [human,mouse]', default = 'human'),
  make_option(c('-o', '--output'), type = 'character', help = 'output', default = 'cell_cycle_results')
)

opt <- parse_args(OptionParser(option_list = option_list))
MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}


library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)