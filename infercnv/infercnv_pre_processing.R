library(optparse)

option_list <- list(
  make_option(c('-s', '--seurat_Obj'), type = 'character', help = 'seurat Object'),
  make_option(c('-g','--group_config'), type = 'character', help = 'set group'),
  make_option(c('-c','--cell_config'), type = 'character', help = 'set ob cells and ref cells')
)
opt <- parse_args(OptionParser(option_list = option_list))


library(tidyverse)
MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}

MyPP <- function(seob, group_config, cell_config){
  MyMkdir('step1_infercnv_pp')
  
  group <- read.table(group_config, sep = '\t', header = F) %>% separate_rows(V2, sep = '@')
  group_anno <- group$V1
  names(group_anno) <- group$V2
  
  MyAddGroup <- function(x){
    if(!is.na(group_anno[x])){
      return(group_anno[[x]]) 
    }else{
      print(str_c(x, 'not in the group file'))
    }
  }
  
  # 读取cell config文件，第二列不可重复
  ref_cell <- read.table(cell_config, sep = '\t', header = F)
  ref_cell_name <- ref_cell$V2
  
  #判断配置文件第二列是否重复
  if(length(unique(ref_cell_name)) != length(ref_cell_name)){
    stop('wrong cell config file!!!')
  }
  
  #判断配置文件第三列是否仅且必须包括ref和ob
  if(length(unique(ref_cell$V3)) != 2){
    stop('wrong cell config file!!!')
  }else if(!('ob' %in% unique(ref_cell$V3) & 'ref' %in% unique(ref_cell$V3))){
    stop('wrong cell config file!!!')
  }
  
  #判断配置文件中是否一个细胞既做了ref也做了ob
  if(dim(ref_cell)[1] != dim(ref_cell %>% distinct(V2, V3))[1]){
    stop('wrong cell config file!!!')
  }

  MyAddCellInfo <- function(sample, cell){
    # 添加标签，
    if(cell %in% ref_cell_name){
      samples <- ref_cell %>% filter(V2 == cell) %>% pull(V1)
      class <- ref_cell %>% filter(V2 == cell) %>% pull(V3)
      if(class == 'ob'){
        class = '_observation'
      }else if(class == 'ref'){
        class = '_reference'
      }else{
        stop('wrong cell config file!!!')
      }
      if(samples == 'all' | sample %in% str_split(samples, '@')[[1]]){
        return(str_c(cell, class))
      }else{
        return(cell)
        }
    }else{
        return(cell)
      }
  }
  
  
  library(Seurat)
  seob <- readRDS(seob)
  
  # meta中添加ref信息，分组信息，降温信息
  print('processing meta')
  seob@meta.data %>% rownames_to_column(var = 'barcodes') %>% 
    rowwise() %>%
    mutate(infercnv_group = MyAddGroup(sample),
           infercnv_reference = MyAddCellInfo(sample, cellType)) %>%
    left_join(as.data.frame(seob@reductions$umap@cell.embeddings) %>%
                rownames_to_column(var = 'barcodes'), 
              by = 'barcodes')  %>%
    filter(str_detect(infercnv_reference, '_observation|_reference')) %>% 
    drop_na() -> meta
  save(meta, file = 'step1_infercnv_pp/meta.RData')
  
  # gene matrix
  print('saving gene matrix')
  counts_matrix <- as.data.frame(seob@assays$RNA@counts)[, meta$barcodes]
  write.table(counts_matrix, file = 'step1_infercnv_pp/raw_counts_matrix.txt', 
              sep = '\t', row.names = T, col.names=T, quote = F);rm(counts_matrix)
  
  # cell_anno
  print('saving cell anno')
  meta %>% select(barcodes, infercnv_reference) %>%
    column_to_rownames(var = 'barcodes') %>%
    write.table(file = 'step1_infercnv_pp/anno_file.txt', sep = '\t', row.names = T, col.names=F, quote = F)
  
  # reference cell name
  print('saving reference cell name')
  meta %>% filter(str_detect(infercnv_reference, '_reference')) %>% 
    distinct(infercnv_reference) %>% pull(infercnv_reference) -> ref_cell
  save(ref_cell, file = 'step1_infercnv_pp/ref_cell.RData')
  
  
  myinfercnv_obj <- list(meta = meta, 
                         gene_matrix = 'step1_infercnv_pp/raw_counts_matrix.txt',
                         cell_anno = 'step1_infercnv_pp/anno_file.txt',
                         ref_cell = ref_cell)
  save(myinfercnv_obj, file = 'step1_infercnv_pp/myinfer_obj.RData')
}


MyPP(seob = opt$seurat_Obj, group_config = opt$group_config, cell_config = opt$cell_config)


