#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c('-s', '--seurat_Obj'), type = 'character', help = 'seurat Object'),
  make_option(c('-p','--species'), type = 'character', 
              help = 'species [human,mouse,rat]', default = 'human'),
  make_option(c('-c','--cell_type'), type = 'character',
              help = 'colnames of seurat Object of cellType [default %default]', default = 'cellType'),
  make_option('--MyMakeMatrix', action = 'store_true', help = 'run MyMakeMatrix'),
  make_option('--MyCpdbPlot', action = 'store_true', help = 'run MyCpdbPlot'),
  make_option(c('-n', '--count_network'), type = 'character', help = 'count_network.txt'),
  make_option(c('-g', '--group_config'), type = 'character', help = 'group config [default %default]', default = 'None'),
  make_option(c('-o', '--output'), type = 'character', help = 'plot output'),
  make_option(c('-r', '--results'), type = 'character', help = 'results [default %default]', default = 'cpdb_resluts')
)

opt <- parse_args(OptionParser(option_list = option_list))
library(tidyverse)

MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}

MyMakeMatrix <- function(SeuratObj,species='human',cell_type='cellType',output='results'){
  "获取seurat中的count矩阵以及matedata"
  MyMkdir(output)

  #读取同源文件
  anno_df <- '/share/nas1/zhangjm/workspace/MyUtils/cpdb/mouse_human_rat_Homologous.txt'
  anno_df <- read.table(anno_df, sep = '\t',header = T) %>% 
    mutate(mouse_gene_symbol = str_to_upper(mouse_gene_symbol),
           rat_gene_symbol = str_to_upper(rat_gene_symbol)) %>%
    select(human_gene_id, all_of(str_c(species, '_gene_symbol'))) %>%
    drop_na()
  
  anno_ve <- anno_df[,1];names(anno_ve) <- anno_df[,2]
  changeID2symbol <- function(gene_id){
    if(!is.na(anno_ve[gene_id])){
      return(anno_ve[[gene_id]])
    }else{
      print(str_c(gene_id, 'is not in the anno_df file'))
      return(NA)
  }}
    
  #导出meta_data
  SeuratObj@meta.data %>% as.data.frame() %>%
    rownames_to_column(var = 'barcodes') %>%
    select(barcodes, all_of(cell_type)) %>%
    drop_na() %>%
    write.table(file = str_c(output, '/meta_data.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
  
  #导出matrix
  gene_count <- SeuratObj@assays$RNA@counts %>% as.data.frame() %>%
    rownames_to_column(var = 'gene_symbol') %>%
    mutate(gene_symbol = str_to_upper(gene_symbol)) %>% rowwise() %>%
    mutate(gene_symbol = changeID2symbol(gene_symbol)) %>%
    drop_na()
  
  # if(species == 'human'){
  #   print('doing human')
  #   gene_count <- SeuratObj@assays$RNA@counts %>% as.data.frame() %>%
  #     rownames_to_column(var = 'gene_symbol') %>%
  #     mutate(gene_symbol = str_to_upper(gene_symbol)) %>%
  #     left_join(anno_df %>% select(human_gene_symbol, human_gene_id),
  #               by = c('gene_symbol'='human_gene_symbol')) %>%
  #     mutate(gene_symbol = human_gene_id) %>%
  #     select(-human_gene_id) %>%
  #     drop_na()
  #   
  # }else if(species == 'mouse'){
  #   print('doing mouse')
  #   gene_count <- SeuratObj@assays$RNA@counts %>% as.data.frame() %>%
  #     rownames_to_column(var = 'gene_symbol') %>%
  #     mutate(gene_symbol = str_to_upper(gene_symbol)) %>%
  #     left_join(anno_df %>% select(mouse_gene_symbol, human_gene_id),
  #               by = c('gene_symbol'='mouse_gene_symbol')) %>%
  #     mutate(gene_symbol = human_gene_id) %>%
  #     select(-human_gene_id) %>%
  #     drop_na()
  #   
  # }else if(species == 'rat'){
  #   print('doing rat')
  #   gene_count <- SeuratObj@assays$RNA@counts %>% as.data.frame() %>%
  #     rownames_to_column(var = 'gene_symbol') %>%
  #     mutate(gene_symbol = str_to_upper(gene_symbol)) %>%
  #     left_join(anno_df %>% select(rat_gene_symbol, human_gene_id),
  #               by = c('gene_symbol'='rat_gene_symbol')) %>%
  #     mutate(gene_symbol = human_gene_id) %>%
  #     select(-human_gene_id) %>%
  #     drop_na()
  # }
  # else{
  #   stop('wrong species! human, mouse, rat are supported.')
  # }
  
  #导出gene matrix文件
  print('writing gene matrix file ...')
  gene_count %>% write.table(file = str_c(output, '/gene_count.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
}


MyCpdbPlot <- function(config, count_network, output){
  "基于igraph绘制网络图"
  count_network <- read_delim(count_network, 
                              delim = "\t", 
                              escape_double = FALSE,
                              trim_ws = TRUE) %>%
    filter(count > 0) %>% rowwise() %>%
    mutate(col = config[SOURCE])
  
  nodes <- data.frame(cell = unique(c(count_network$SOURCE,
                                      count_network$TARGET))) %>%
    mutate(col = config[cell])
  
  net<- graph_from_data_frame(d = count_network, vertices = nodes, directed = T)
  E(net)$width  <- E(net)$count/50
  V(net)$color <- nodes$col
  E(net)$color <- count_network$col
  coords <- layout_in_circle(net, order = names(config))
  
  #修改颜色
  # all_color = getPalette(celltype_num)
  # for (i in 1: length(unique(count_network$SOURCE)) ){
  #   E(net)[map(unique(count_network$SOURCE),function(x) {
  #     get.edge.ids(net,vp = c(unique(count_network$SOURCE)[i],x))
  #   })%>% unlist()]$color <- all_color[i]
  # }
  
  png(str_c(output, '/cell_talk.png'), res = 1200, width = 5000, height = 5000)
  plot(net, edge.arrow.size=.2, 
       edge.curved= .2,
       #vertex.color=getPalette(celltype_num),
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=.45, width=0.1) 
  dev.off()
  
  pdf(str_c(output, '/cell_talk.pdf'))
  plot(net, edge.arrow.size=.2, 
       edge.curved= .2,
       #vertex.color=getPalette(celltype_num),
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=.45, width=0.1)
  dev.off()
}


### 准备cellphonedb的输入数据
if(isTRUE(opt$MyMakeMatrix)){
  library(Seurat)
  
  if(!opt$species %in% c('human','mouse','rat')){
    stop('wrong species! human, mouse, rat are supported.')
  }
  
  if(opt$group_config == 'None'){
    SeuratObj <- readRDS(opt$seurat_Obj)
    des = str_c(opt$results, '/', 'all_sample_cpdb')
    MyMakeMatrix(SeuratObj = SeuratObj, 
                 species = opt$species, 
                 cell_type = opt$cell_type, 
                 output = des)
    
  }else{
    # 解析分组config信息
    group <- read.table(opt$group_config, sep = '\t', header = F) %>% separate_rows(V2, sep = '@')
    group_anno <- group$V1;names(group_anno) <- group$V2
    
    #添加分组
    MyAddGroup <- function(x){
      if(!is.na(group_anno[x])){
        return(group_anno[[x]]) 
      }else{
        print(str_c(x, ' not in the group config file'))
        return('NA')
      }
    }
    
    SeuratObj <- readRDS(opt$seurat_Obj)
    SeuratObj@meta.data %>% rownames_to_column(var = 'barcodes') %>%
      rowwise() %>%
      mutate(cpdb_group = MyAddGroup(sample)) %>% 
      column_to_rownames(var = 'barcodes') -> SeuratObj@meta.data
    
    
    # 按照指定列拆分
    SeuratObj <- subset(SeuratObj, subset = cpdb_group != 'NA')
    seob_list <- SplitObject(SeuratObj, split.by='cpdb_group')
    
    for(file_name in names(seob_list)){
      des = str_c(opt$results, '/', file_name, '_cpdb')
      MyMakeMatrix(SeuratObj = seob_list[[file_name]], 
                   species = opt$species, 
                   cell_type = opt$cell_type, 
                   output = des)
    }
  }
}

### 基于cellphonedb结果绘图
if(isTRUE(opt$MyCpdbPlot)){
  library(RColorBrewer)
  library(tidyverse)
  library(igraph)
  getPalette = colorRampPalette(brewer.pal(8, "Set2"))
  
  count_network <- read_delim(opt$count_network, 
                              delim = "\t", 
                              escape_double = FALSE,
                              trim_ws = TRUE)
  
  cells <- unique(c(count_network$SOURCE, count_network$TARGET)) %>% sort()
  col_ve = getPalette(length(cells))
  names(col_ve) <- cells
  MyCpdbPlot(config = col_ve,count_network = opt$count_network, output = opt$output)
  }
