.libPaths('/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/lib/R/library') 
library(optparse)

option_list <- list(
  make_option(c('-s', '--seurat_Obj'), type = 'character', help = 'seurat Object'), 
  make_option(c('-o', '--output'), type = 'character', help = 'plot output', default = 'pySCENIC_results')
)

opt <- parse_args(OptionParser(option_list = option_list))

get_loom_and_metadata <- function(seob_obj, output){
  library(SeuratDisk)
  library(Seurat)
  library(tidyverse)
  MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}
  MyMkdir(output)
  
  seob <- readRDS(seob_obj)
  meta_data <- seob@meta.data %>% rownames_to_column('barcodes')
  save(meta_data, file = file.path(output, 'meta_data.RData'))
  
  DefaultAssay(seob) <- 'RNA'
  my_loom <- as.loom(seob, filename=file.path(output, 'seob_obj.loom'))
  my_loom$close_all()
  
}

get_loom_and_metadata(seob_obj = opt$seurat_Obj, output = opt$output)
