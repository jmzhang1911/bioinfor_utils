.libPaths('/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/lib/R/library') 
library(optparse)

option_list <- list(
  make_option(c('-s', '--seob_obj'), type = 'character', help = 'Seurat Obj'), 
  make_option(c('-o', '--output'), type = 'character', help = 'plot output', default = 'GSVA_results'),
  make_option(c('-c', '--cell_type'), type = 'character', help = 'colnames of celltype', default = 'cellType'),
  make_option(c('-p', '--species'), type = 'character', help = 'species [hs|mm]', default = 'hs'),
  make_option(c('-g', '--gs_cat'), type = 'character', help = '[C1|C2|C3|C4|C5|C6|C7|C8|H], sep by ,', default = 'all')
)

opt <- parse_args(OptionParser(option_list = option_list))


suppressMessages({
  library(funr)
  library(msigdbr)
  library(Seurat)
  library(GSVA)
  library(tidyverse)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(patchwork)
  library(RColorBrewer)
  library(circlize)
  library(cowplot)
  library(ggplotify)
  library(foreach)
  library(doParallel)
  library(ggsci)
  
  source(file.path(dirname(sys.script()), '../Utils.R'))
})


MyGSVAPlot <- function(matrix, marker_genes, title = 'GSVA results Top5'){
  Do_col_aut <- function(){
    # find the default legend scale and color it
    heatmap_col = c('#6BAED6', 'white','#FDBB84')
    p_res <- Heatmap(matrix)
    max_col <- max(p_res@matrix_color_mapping@levels)
    min_col <- min(p_res@matrix_color_mapping@levels)
    col = colorRamp2(breaks = c(min_col, mean(c(max_col, min_col)), max_col),
                     colors = heatmap_col)
    return(col)
  }
  
  row_anno <- rowAnnotation(foo = anno_mark(at = match(marker_genes, rownames(matrix[marker_genes, ])), 
                                            labels = str_replace_all(marker_genes, '_', ' ') %>%
                                              str_wrap(width = 60), 
                                            labels_gp = gpar(fontsize = 5)))
  
  Heatmap(matrix[marker_genes, ], col = Do_col_aut(),
          show_row_dend = T, show_column_dend = T,
          column_dend_height = unit(.3, 'cm'),
          right_annotation = row_anno,
          show_row_names = F,
          show_column_names = T, 
          column_title = title,
          column_names_gp = gpar(fontsize = unit(8, 'npc')),
          row_names_gp = gpar(fontsize = unit(6, 'npc')),
          raster_quality = 30, 
          heatmap_legend_param = list(title = 'score',
                                      title_gp = gpar(fontsize = unit(8, 'npc'),
                                                      fontface = 'bold'),
                                      legend_height = unit(0.3, 'npc'),
                                      title_position = 'leftcenter-rot'),
          width = 18, height = 12) -> p
  
  Heatmap(matrix, col = Do_col_aut(),
          show_row_dend = T, show_column_dend = T,
          column_dend_height = unit(.3, 'cm'),
          right_annotation = row_anno,
          show_column_names = T, 
          show_row_names = F,
          column_title = title,
          column_names_gp = gpar(fontsize = unit(8, 'npc')),
          row_names_gp = gpar(fontsize = unit(6, 'npc')),
          raster_quality = 30, 
          heatmap_legend_param = list(title = 'score',
                                      title_gp = gpar(fontsize = unit(8, 'npc'),
                                                      fontface = 'bold'),
                                      legend_height = unit(0.3, 'npc'),
                                      title_position = 'leftcenter-rot'),
          width = 18, height = 12) -> p2
  
  
  return(list(p=p, p2=p2))
}


run <- function(term, gene_expr, output = './'){
  # 结果路径
  term %>%
    distinct(results_path) %>%
    pull(results_path) -> terms_path
  results_path <- file.path(output, terms_path)
  MyMkdir(results_path)
  
  # 查分term
  split.data.frame(term, term$gs_name) %>% 
    lapply(function(x){return(x$gene_symbol)}) -> tmp_gene_set
  
  # 运算GSVA
  gsva_results <- tryCatch({
    gsva(as.matrix(gene_expr), tmp_gene_set)
  }, error = function(e) {
  })
  
  if(is.null(gsva_results)){
    return()
  }
  
  nrow_number <- nrow(as.data.frame(gsva_results))
  if( nrow_number< 1){
    print('no results')
    return(NULL)
  }
  
  gsva_results %>% as.data.frame() %>% 
    rownames_to_column('terms') %>%
    pivot_longer(cols = -1, names_to = 'cellType', values_to = 'value') -> all_terms_results 
  
  # 去TOP5
  all_terms_results %>% 
    group_by(cellType) %>% 
    slice_max(value, n = 5) %>% 
    pull(terms) %>% unique() -> top5_terms
  
  
  # 保存热图
  
  title <- str_c('category of ', terms_path) %>% str_replace_all('\\/', ' ') %>% str_remove_all('\\.')
  MyGSVAPlot(gsva_results, top5_terms, title = title) -> plot_list
  
  pdf(file = file.path(results_path, 'all_terms_heatmap.pdf'))
  draw(plot_list$p2,  heatmap_legend_side = "left")
  dev.off()
  
  png(file = file.path(results_path, 'all_terms_heatmap.png'), width = 1200, height = 1200, res = 200)
  draw(plot_list$p2,  heatmap_legend_side = "left")
  dev.off()
  
  pdf(file = file.path(results_path, 'top5_terms_heatmap.pdf'))
  draw(plot_list$p,  heatmap_legend_side = "left")
  dev.off()
  
  png(file = file.path(results_path, 'top5_terms_heatmap.png'), width = 1200, height = 1200, res = 200)
  draw(plot_list$p,  heatmap_legend_side = "left")
  dev.off()
  
  # 保存所有结果
  write.table(all_terms_results, 
              file = file.path(results_path, 'all_terms_results.txt'), 
              sep = '\t', quote = F, row.names = F)
}

make_summary(sys.script(), 'doing')
if(T){
  
  if(opt$species == 'hs'){
    load(file.path(dirname(sys.script()), 'hs_msigdbr_data.RData'))
  }else if(opt$species == 'mm'){
    load(file.path(dirname(sys.script()), 'mm_msigdbr_data.RData'))
  }else{
    stop('wrong species, only for hs|mm') 
  }
  
  print('making dataset ...')
  df <- msigdbr_data %>% ungroup()
  if (opt$gs_cat != 'all') {
    gs_cat_vector <- str_split(opt$gs_cat, ',')[[1]]
    df <- df %>% filter(gs_cat %in% gs_cat_vector)
    print(str_c('using ', gs_cat_vector))
  }
  
  split.data.frame(df, df$results_path) -> all_gene_set_list
  
  print('getting expr matrix')
  readRDS(opt$seob_obj) %>% 
    AverageExpression(assays = 'RNA', slot = 'data', group.by = opt$cell_type) -> expr
  expr <- expr$RNA[rowSums(expr$RNA) > 0, ]
  
  MyMkdir(opt$output)
  write.table(expr, file = file.path(opt$output, 'gene_expr.txt'), sep = '\t', row.names = T, quote = F)
  
  gc()
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  print('doing GSVA ...')
  foreach(x = all_gene_set_list,
          .packages = c('tidyverse', 'enrichplot', 'GSVA', 'ComplexHeatmap',
                        'ggplotify', 'clusterProfiler', 'circlize',
                        'parallel', 'doParallel'),
          .export = c('MyMkdir', 'MyGSVAPlot')) %dopar%
    run(term = x, gene_expr = expr, output = opt$output)
  stopCluster(cl)
}
make_summary(sys.script(), 'done')

####测试
# read_delim("gene_expr.txt", 
#            delim = "\t", escape_double = FALSE, 
#            trim_ws = TRUE) %>% column_to_rownames('xxx') -> gene_expr
# load('msigdbr_data.RData') # 加载整理后的基因集
# split.data.frame(msigdbr_data, msigdbr_data$results_path) -> all_gene_set_list
# 
# term <- all_gene_set_list$`C7/VAX`
# split.data.frame(term, term$gs_name) %>% 
#   lapply(function(x){return(x$gene_symbol)}) -> tmp_gene_set
# gsva_results <- gsva(as.matrix(gene_expr), tmp_gene_set) 
# gsva_results %>% as.data.frame() %>% 
#   rownames_to_column('terms') %>%
#   pivot_longer(cols = -1, names_to = 'cellType', values_to = 'value')  %>% 
#   group_by(cellType) %>% 
#   slice_max(value, n = 5) %>% 
#   pull(terms) %>% unique() -> top5_terms
# 
# 
# matrix <- as.matrix(gsva_results)
# 
# title = 'sb'
# Heatmap(matrix, 
#         show_row_dend = T, 
#         show_column_dend = T,
#         column_dend_height = unit(.3, 'cm'),
#         #row_names_max_width = unit(5, "cm"),
#         #row_labels = str_wrap(rownames(matrix), width = 50),
#         show_column_names = T, 
#         right_annotation = row_anno,
#         show_row_names = F,
#         column_title = title,
#         column_names_gp = gpar(fontsize = unit(8, 'npc')),
#         row_names_gp = gpar(fontsize = unit(6, 'npc')),
#         raster_quality = 30, 
#         heatmap_legend_param = list(title = 'score',
#                                     title_gp = gpar(fontsize = unit(8, 'npc'),
#                                                     fontface = 'bold'),
#                                     legend_height = unit(0.3, 'npc'),
#                                     title_position = 'leftcenter-rot'),
#         width = 12, height = 12) -> p
# 
# 
# run(gene_expr = gene_expr, term = term)
# 
# term %>%
#   distinct(results_path) %>%
#   pull(results_path) -> results_path

