
library(optparse)

option_list <- list(
  make_option(c('-s', '--seurat_Obj'), type = 'character', help = 'seurat Object'),
  make_option(c('-p','--species'), type = 'character', 
              help = 'species [human,mouse]', default = 'human'),
  make_option(c('-o', '--output'), type = 'character', help = 'output', default = 'cell_cycle_results')
)

opt <- parse_args(OptionParser(option_list = option_list))
MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}


library(tidyverse)
library(Seurat)


cell_cycle_score <- function(seob, species, output){
  suppressPackageStartupMessages(library(here))
  #here::i_am('./run_cellcycle.R');pwd <- here()
  pwd <- '/share/nas1/zhangjm/workspace/MyUtils/cell_cycle_analysis/'
  
  # 拷贝readme文件
  file.copy(file.path(pwd, './cell_cycle_readme.zip'), to = output)
  
  MyMkdir(output)
  
  print('reading rds ...')
  seob <- readRDS(seob)
  seob_cycle <- CreateSeuratObject(counts = seob@assays$RNA@counts, meta.data = seob@meta.data)
  
  seob_cycle %<>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst") %>%
    ScaleData() 
  

  if(species == 'human'){
    load(str_c(pwd, '/Human_cell_cycle_gene.rda'))
  }else if(species == 'mouse'){
    load(str_c(pwd, '/Mouse_cell_cycle_gene.Rda'))
  }else{
    stop('only support mouse and human')
  }
  
  seob_cycle <- CellCycleScoring(seob_cycle, s.features = s_genes, g2m.features = g2m_genes)

  seob_cycle@meta.data %>% rownames_to_column(var = 'barcodes') %>% 
    select(barcodes, S.Score, G2M.Score, Phase) -> cycle_df;rm(seob_cycle)
  
  seob@meta.data %<>% rownames_to_column(var = 'barcodes') %>%
    left_join(cycle_df, by = 'barcodes') %>% 
    column_to_rownames(var = 'barcodes')
  
  write.table(seob@meta.data, file = str_c(output, '/cell_cycle.txt'), col.names = T, row.names = T, sep = '\t', quote = F)
  print('saving rds ...')
  saveRDS(seob, file = str_c(output, '/seob_cell_cycled.Rds'))
  res <- list(seob = seob, s_genes = s_genes, g2m_genes = g2m_genes)
  return(res)
}


sb_plot <- function(seob, od){
  library(pheatmap)
  print('doing sb plotting ...')
  prefix <- 'seob_cell_cycled'
  
  seurat.obj <- seob$seob
  s_genes <- seob$s_genes
  g2m_genes <- seob$g2m_genes
  
  p <- DimPlot(seurat.obj, reduction = "umap", group.by = "Phase")
  P.umap <- UMAPPlot(seurat.obj, group.by = "Phase", label = FALSE) + labs(title = "")
  P.tsne <- TSNEPlot(seurat.obj, group.by = "Phase", label = FALSE) + labs(title = "")
  ggsave(file = file.path(od, paste(prefix, "cellCycle_umap.png", sep = ".")), P.umap, width = 6, height = 4, scale = 1.3)
  ggsave(file = file.path(od, paste(prefix, "cellCycle_umap.pdf", sep = ".")), P.umap, width = 6, height = 4, scale = 1.3)
  ggsave(file = file.path(od, paste(prefix, "cellCycle_tsne.png", sep = ".")), P.tsne, width = 6, height = 4, scale = 1.3)
  ggsave(file = file.path(od, paste(prefix, "cellCycle_tsne.pdf", sep = ".")), P.tsne, width = 6, height = 4, scale = 1.3)
  
  CC <- seurat.obj@meta.data %>% rownames_to_column(var = "Cell") %>% select(Cell, sample, S.Score, G2M.Score, Phase)
  write.table(CC, file.path(od, paste(prefix, "cell_cycle.xls", sep = ".")), col.names = T, row.names = F, quote = F, sep = "\t")
  sample_id <- unique(CC$sample)
  CC_cluster <- seurat.obj@meta.data %>% rownames_to_column(var = "Cell") %>% select(Cell, sample, Phase, seurat_clusters)
  
  ###########pheatmap
  ##col_ann
  col_ann <- CC %>% select(Cell, Phase)
  col_ann <- col_ann[order(col_ann$Phase),] 
  rownames(col_ann) <- col_ann$Cell
  col_ann <- col_ann %>% select(-Cell)
  
  ###data
  all_data <- data.frame(seurat.obj@assays$RNA@data, check.names = F)
  s_p <- all_data %>% rownames_to_column(var = "gene") %>% filter(gene %in% c(s_genes, g2m_genes)) %>% column_to_rownames(var = "gene")
  ##row_ann
  row_ann <- data.frame(gene=rownames(s_p), Class = "")
  row_ann <- row_ann %>% mutate(Class = case_when(gene %in% g2m_genes~"G2/M", gene %in% s_genes~"S")) 
  row_ann <- row_ann[order(row_ann$Class),]
  rownames(row_ann) <- row_ann$gene
  row_ann <- row_ann %>% select(-gene)
  
  s_p <- s_p[,rownames(col_ann)]
  s_p <- s_p[rownames(row_ann),]
  
  ##设置颜色条
  #breaks
  bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
  
  pheatmap(s_p, scale = "row", color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
           legend_breaks=seq(-2, 2, 1), breaks = bk,
           show_rownames = F, show_colnames = F,
           cluster_cols = F, cluster_rows = F,
           annotation_col = col_ann,angle_col = "90",filename = file.path(od, paste(prefix, "cell_cycle_heatmap.pdf", sep = ".")))
  pheatmap(s_p, scale = "row", color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
           legend_breaks=seq(-2, 2, 1), breaks = bk,
           show_rownames = F, show_colnames = F,
           cluster_cols = F, cluster_rows = F,
           annotation_col = col_ann, angle_col = "90", filename = file.path(od, paste(prefix, "cell_cycle_heatmap.png", sep = ".")))
  }


seob <- cell_cycle_score(opt$seurat_Obj, opt$species, opt$output)
sb_plot(seob, opt$output)


