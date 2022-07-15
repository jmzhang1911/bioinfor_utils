.libPaths('/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/lib/R/library') 
library(optparse)

option_list <- list(
  make_option(c('-m', '--meta_data'), type = 'character', help = 'meta_data.RData'), 
  make_option(c('-o', '--output'), type = 'character', help = 'dir name of outputs', default = 'pySCENIC_results'),
  make_option(c('-c', '--cell_type'), type = 'character', help = 'colnames of celltype', default = 'cellType'),
  make_option(c('-g', '--groups'), type = 'character', 
              help = 'colnames to anno in heatmap', 
              default = 'orig.ident'),
  make_option(c('-l', '--sample_loom'), type = 'character', help = 'loom results from pySCENIC'))

opt <- parse_args(OptionParser(option_list = option_list))

suppressMessages({
  library(SCENIC)
  library(tidyverse)
  library(SCopeLoomR)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(pheatmap)
  library(circlize)
  library(Cairo)
})


MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}}

MyPySenicHeatPlot <- function(expr, meta_data, cell_type='cellType', groups='orig.ident', output='pySCENIC_results'){
  # dplyr无法直接从函数参数重获取
  cell_type = cell_type
  
  getPalette_cell = colorRampPalette(brewer.pal(8, "Set2"))
  getPalette_sample = colorRampPalette(brewer.pal(8, "Accent"))
  
  #设置注释颜色
  meta_data %>% 
    filter(barcodes %in% colnames(expr)) %>%
    drop_na() -> meta_df
  
  getPalette_cell(length(unique(meta_df[[cell_type]]))) -> cells
  meta_df[[cell_type]] %>% unique() -> names(cells)
  
  getPalette_sample(length(unique(meta_df[[groups]]))) -> samples
  meta_df[[groups]] %>% unique() -> names(samples)
  
  print(cells);print(samples)
  
  #获取每个细胞类型中最活跃的转录因子
  expr %>% t() %>% as.data.frame() %>% rownames_to_column('barcodes') %>% 
    filter(barcodes %in% meta_df$barcodes) %>%
    pivot_longer(cols = -1, names_to = 'TF', values_to = 'value') %>%
    left_join(meta_df %>% select(barcodes, !!as.symbol(cell_type)), by = 'barcodes') %>%
    group_by(!!as.symbol(cell_type), TF) %>%
    summarise(mean_value = mean(value)) %>% 
    slice_max(mean_value, n = 5) %>% pull(TF) %>% unique() -> marker
  
  print(marker)
  
  # 极高和极低值处理
  expr <- expr[,colnames(expr) %in% meta_df$barcodes]
  expr[expr > 2] <- 2 
  expr[expr < -2] <- -2
  

  
  
  # 列注释
  top_anno_ref <- HeatmapAnnotation(
    samples = meta_df[[groups]], 
    celltype =  meta_df[[cell_type]], 
    height = unit(0.1, 'npc'),
    col = list(celltype = cells, samples = samples),
    show_legend = FALSE
  )
  
  
  # 行注释
  # 调整大小
  if(length(marker) >= 50){
    fontsize <- 3
  }else{
    fontsize <- 6
  }
  
  row_anno <- rowAnnotation(foo = anno_mark(at = match(marker, rownames(expr)), 
                                            labels = marker, 
                                            labels_gp = gpar(fontsize = fontsize)))
  row_anno2 <- rowAnnotation(foo = anno_mark(at = match(marker, rownames(expr[marker, ])), 
                                            labels = marker, 
                                            labels_gp = gpar(fontsize = fontsize)))
  
  # 整体颜色
  col_fun = colorRamp2(c(-2, 0, 2), c("#87CEFA", "white", "#CC2121")) 
  
  # 主题
  Heatmap(expr, 
          row_names_gp = gpar(fontsize = 1),
          col = col_fun, 
          border = TRUE,
          row_gap = unit(0, "mm"), 
          column_gap = unit(0, "mm"),
          column_title = NULL,
          column_split  = meta_df[cell_type],
          column_title_gp  = gpar(fontsize = 8), 
          top_annotation = top_anno_ref, 
          right_annotation = row_anno,
          show_column_names = F,
          show_row_names = F,
          show_row_dend = F,
          show_column_dend = F,
          show_heatmap_legend = FALSE,
          raster_quality = 15,
          width = 10, height = 14) ->  p
  
  print(dim(expr[marker,]))
  
  
  
  Heatmap(expr[marker,], 
          row_names_gp = gpar(fontsize = 1),
          col = col_fun, 
          border = TRUE,
          row_gap = unit(0, "mm"), 
          column_gap = unit(0, "mm"),
          column_title = NULL,
          column_split  = meta_df[cell_type],
          column_title_gp  = gpar(fontsize = 8), 
          top_annotation = top_anno_ref, 
          right_annotation = row_anno2,
          show_column_names = F,
          show_row_names = F,
          show_row_dend = F,
          show_column_dend = F,
          show_heatmap_legend = FALSE,
          raster_quality = 15,
          width = 10, height = 14) ->  p1
  
  # 图例处理
  legend_expr <- Legend(col_fun = col_fun, direction = 'horizontal')
  legend_cells <-  Legend(labels = names(cells),
                          legend_gp = gpar(fill = as.character(cells)), 
                          title = 'cellType', 
                          direction = "horizontal")
  legend_samples <- Legend(labels = names(samples), 
                           title = 'samples',
                           legend_gp = gpar(fill = as.character(samples)))
  pd <-  packLegend(legend_cells,
                    legend_expr,
                    legend_samples)
  
  MyMkdir(output)
  
  pdf(file.path(output, 'SCENIC_analysis.pdf'))
  draw(p, annotation_legend_list = pd, column_title = 'SCENIC analysis')
  dev.off()
  
  png(file.path(output, 'SCENIC_analysis.png'), height = 1200, width = 1100, res = 200)
  draw(p, annotation_legend_list = pd, column_title = 'SCENIC analysis')
  dev.off()
  
  
  pdf(file.path(output, 'SCENIC_analysis_top5.pdf'))
  draw(p1, annotation_legend_list = pd, column_title = 'SCENIC analysis')
  dev.off()
  
  png(file.path(output, 'SCENIC_analysis_top5.png'), height = 1200, width = 1100, res = 200)
  draw(p1, annotation_legend_list = pd, column_title = 'SCENIC analysis')
  dev.off()
  
}


if(T){
  print('loading data ...')
  
  load(opt$meta_data)
  loom <- open_loom(opt$sample_loom)
  
  
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC")
  raw_matrix <- t(scale(t(getAUC(regulonAUC[,]))))
  
  print('plotting heatmap ...')
  MyPySenicHeatPlot(expr = raw_matrix,
                    meta_data = meta_data, 
                    groups = opt$groups, 
                    cell_type = opt$cell_type, 
                    output = opt$output)
}


