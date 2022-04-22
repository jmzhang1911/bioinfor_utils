library(optparse)
library(tidyverse)

option_list <- list(
  make_option(c('-i', '--myinfercnv_obj'), type = 'character', help = 'myinfercnv_obj'),
  make_option(c('-t','--thresholds'), type = 'character', 
              help = 'thresholds, ps 0.8,1,1.2 [default auto]',
              default = 'auto')
)
opt <- parse_args(OptionParser(option_list = option_list))

MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggsci)
library(Seurat)
getPalette_chr = colorRampPalette(brewer.pal(12, "Set3"))
getPalette_ref = colorRampPalette(brewer.pal(8, "Set2"))
getPalette_ob = colorRampPalette(brewer.pal(9, "Set1"))
getPalette_sample = colorRampPalette(brewer.pal(8, "Accent"))

MyHeatMapPlot <- function(myinfercnv_obj, output){
  infer_obj <- myinfercnv_obj$infer_obj
  subclusters_cell_groupings <- myinfercnv_obj$subclusters_cell_groupings
  meta_df <- myinfercnv_obj
  group_colname <- 'infercnv_group'
  

  #原始对象
  raw_data <- readRDS(infer_obj)
  #参考细胞
  ref_cell_names <- names(raw_data@reference_grouped_cell_indices)
  GetName <- function(x){
    tmp <- str_split(x, '\\.')[[1]]
    if(tmp[1] %in% ref_cell_names){
      return(tmp[1])
    }else{
      return(str_c(tmp[2:length(tmp)], collapse = '.'))
    }
  }
  
  #所有细胞及亚克隆
  sub_barcodes <- read.table(subclusters_cell_groupings, sep = '\t', header = T) %>% 
    rowwise() %>% mutate(subcluster = GetName(cell_group_name))
  
  # 原始的所有barcodes，gene矩阵
  gene_matrix <- raw_data@expr.data %>% as.data.frame()
  
  ##################### 观测细胞热图
  # 热图矩阵
  matrix_ob <- gene_matrix[,filter(sub_barcodes, !subcluster %in% ref_cell_names) %>% pull(cell)]
  
  # 热图整体颜色
  #Heatmap(t(gene_matrix), cluster_rows = F, cluster_columns = F)@matrix_color_mapping@levels -> legend_scale
  if(opt$thresholds == 'auto'){
    legend_scale <- read.table(myinfercnv_obj$thresholds, header = F)$V1
    if(round(min(legend_scale), 1) > min(legend_scale)){
      min_legend <- round(min(legend_scale), 1) - 0.1
    }else{
      min_legend <- round(min(legend_scale), 1)
    }
    
    if(round(max(legend_scale), 1) > max(legend_scale)){
      max_legend <- round(max(legend_scale), 1)
    }else{
      max_legend <- round(max(legend_scale), 1) + 0.1
    }  
    
    mid_legend <- round(min_legend + (max_legend - min_legend) / 2 , 1)
    col_fun = colorRamp2(c(min_legend, mid_legend, max_legend), c("#000080", "white", "#B22222"))
  }else{
    min_legend = as.numeric(opt$thresholds[1])
    mid_legend = as.numeric(opt$thresholds[2])
    max_legend = as.numeric(opt$thresholds[3])
    col_fun = colorRamp2(c(min_legend, mid_legend, max_legend), c("#000080", "white", "#B22222"))
  }

  
  # 基因信息
  gene_anno_ob <- raw_data@gene_order %>% 
    rownames_to_column(var = 'gene_symbol') %>% 
    filter(gene_symbol %in% rownames(matrix_ob))
  # 细胞信息
  cell_anno_ob <- sub_barcodes %>% filter(!subcluster %in% ref_cell_names) 
  write.table(cell_anno_ob, file = str_c(output, '/subcluster_infor.txt'), sep = '\t', quote = F, row.names = F)
  
  # 生成一个样本信息
  meta <- myinfercnv_obj$meta %>% select(barcodes, group = all_of(group_colname), UMAP_1, UMAP_2)
  cell_anno_ob <- left_join(cell_anno_ob, meta %>% mutate(group = as.character(group)),
                            by = c('cell'='barcodes')) 
  
  top_anno_ob <- HeatmapAnnotation(
    pos = gene_anno_ob$chr, height = unit(0.1, 'npc'),
    show_legend = F, annotation_label = NULL
  )
  
  getPalette_ob(length(unique(cell_anno_ob$subcluster))) -> ob_cells
  cell_anno_ob$subcluster %>% unique() -> names(ob_cells)
  
  getPalette_sample(length(unique(cell_anno_ob$group))) -> ob_group
  cell_anno_ob$group %>% unique() -> names(ob_group)
  
  left_anno_ob <- rowAnnotation(
    `observations cells` = cell_anno_ob$subcluster,
    col = list(`observations cells` = ob_cells),
    show_legend = F, annotation_label = ''
  )
  
  right_anno_ob <- rowAnnotation(
    groups = cell_anno_ob$group,
    col = list(groups = ob_group), 
    show_legend = F
  )
  
  legend_expr <-  Legend(col_fun = col_fun, title = "Modified Expression", direction = "horizontal")
  legend_ob_cells <-  Legend(labels = names(ob_cells) %>% str_remove_all('_observation'), size = 1,
                             legend_gp = gpar(fill = as.character(ob_cells)), 
                             title = "observations cells", ncol = 2, gap = unit(0.86, 'cm')) 
  legend_ob_group <- Legend(labels = names(ob_group), 
                            legend_gp = gpar(fill = as.character(ob_group)), 
                            title = "groups")
  # 绘制观测热图
  Heatmap(t(matrix_ob),
          show_heatmap_legend = F,
          column_dend_reorder = T,
          col = col_fun,
          border = TRUE,
          row_gap = unit(0, "mm"), 
          column_gap = unit(0, "mm"),
          column_split = factor(gene_anno_ob$chr), 
          row_split = factor(cell_anno_ob$group),
          show_column_names = F,
          left_annotation = left_anno_ob, 
          right_annotation = right_anno_ob,
          #top_annotation = top_anno_ob, 
          column_title_rot = 90,
          column_title_gp  = gpar(fontsize = 10),
          row_title = 'observations cells',
          show_row_names = F, 
          cluster_columns = F,
          cluster_rows = F, height = unit(18, 'cm'), width = unit(30, 'cm')) -> p_ob
  
  
  ######### 参考细胞热图
  # 热图矩阵
  matrix_ref <- gene_matrix[,filter(sub_barcodes, subcluster %in% ref_cell_names) %>% pull(cell)]
  # 基因信息
  gene_anno_ref <- raw_data@gene_order %>% 
    rownames_to_column(var = 'gene_symbol') %>% 
    filter(gene_symbol %in% rownames(matrix_ref))
  # 细胞信息
  cell_anno_ref <- sub_barcodes %>% filter(subcluster %in% ref_cell_names)
  
  # 染色体颜色
  chr_ve <- unique(gene_anno_ref$chr)
  chr_ve <- getPalette_chr(length(chr_ve))
  names(chr_ve) <- unique(gene_anno_ref$chr)
  
  top_anno_ref <- HeatmapAnnotation(
    pos = gene_anno_ref$chr,
    height = unit(0.01, 'npc'),
    show_legend = F,
    col = list(pos = chr_ve), annotation_label = ''
  )
  
  getPalette_ref(length(unique(cell_anno_ref$subcluster))) -> ref_cells
  cell_anno_ref$subcluster %>% unique() -> names(ref_cells)
  
  
  left_anno_ref <- rowAnnotation(
    `reference cells` = cell_anno_ref$subcluster,
    annotation_name_side = 'top',
    col = list(`reference cells` = ref_cells),
    show_legend = F, annotation_label = ''
  )
  
  legend_ref_cells <-  Legend(labels = names(ref_cells), 
                              legend_gp = gpar(fill = as.character(ref_cells)), 
                              title = "reference cells")
  
  Heatmap(t(matrix_ref),
          show_heatmap_legend = F,
          column_dend_reorder = T,
          col = col_fun,
          border = TRUE,
          row_gap = unit(0, "mm"), 
          column_gap = unit(0, "mm"),
          column_split = factor(gene_anno_ref$chr), 
          row_split = factor(cell_anno_ref$subcluster),
          show_column_names = F,
          left_annotation = left_anno_ref, 
          top_annotation = top_anno_ref, 
          column_title_rot = 90,
          column_title_gp  = gpar(fontsize = 8),
          row_title = 'reference cells',
          show_row_names = F, 
          cluster_columns = F,
          cluster_rows = F, height = unit(2, 'cm'), width = unit(30, 'cm')) -> p_ref
  
  pd <-  packLegend(legend_expr, 
                    legend_ref_cells,
                    legend_ob_group, 
                    direction = 'horizon')
  
  pd_all <- packLegend(pd, legend_ob_cells)
  
  # 清理以下内存
  rm(raw_data);gc()
  
  ### 添加分组信息但未对subcluster聚类
  png(filename = str_c(output, '/infercnv_subclusters_heatmap.png'), width = 1200, height = 800)
  draw(p_ref %v% p_ob, annotation_legend_side = 'right',
       annotation_legend_list = pd_all, column_title = 'inferCNV subclusters', 
       column_title_gp = gpar(fontsize = 25),
       use_raster = TRUE, raster_quality = 30)
  dev.off()
  
  pdf(str_c(output, '/infercnv_subclusters_heatmap.pdf'), width = 120, height = 80)
  draw(p_ref %v% p_ob, annotation_legend_side = 'right',
       annotation_legend_list = pd_all, column_title = 'inferCNV subclusters', 
       column_title_gp = gpar(fontsize = 25),
       use_raster = TRUE, raster_quality = 30)
  dev.off()
  
  ### 不添加分组信息但对subcluster做聚类
  Heatmap(t(matrix_ob),
          show_heatmap_legend = F,
          column_dend_reorder = T,
          col = col_fun,
          border = TRUE,
          row_gap = unit(0, "mm"), 
          column_gap = unit(0, "mm"),
          column_split = factor(gene_anno_ob$chr), 
          row_split = factor(cell_anno_ob$subcluster),
          show_column_names = F,
          left_annotation = left_anno_ob, 
          column_title_rot = 90,
          column_title_gp  = gpar(fontsize = 10),
          row_title = 'observations cells',
          show_row_names = F, 
          cluster_columns = F,
          cluster_rows = T, 
          height = unit(18, 'cm'), 
          row_dend_width = unit(1.5, 'cm'),
          width = unit(30, 'cm')) -> p_ob_clustered
  
  pd <-  packLegend(legend_expr, 
                    legend_ref_cells,
                    direction = 'horizon')
  
  pd_all <- packLegend(pd, legend_ob_cells)
  
  png(filename = str_c(output, '/infercnv_subclusters_heatmap_clustered.png'), width = 1400, height = 800)
  draw(p_ref %v% p_ob_clustered, annotation_legend_side = 'right',
       annotation_legend_list = pd_all, column_title = 'inferCNV subclusters', 
       column_title_gp = gpar(fontsize = 25),
       use_raster = TRUE, raster_quality = 30)
  dev.off()
  
  pdf(str_c(output, '/infercnv_subclusters_heatmap_clustered.pdf'), width = 120, height = 80)
  draw(p_ref %v% p_ob_clustered, annotation_legend_side = 'right',
       annotation_legend_list = pd_all, column_title = 'inferCNV subclusters', 
       column_title_gp = gpar(fontsize = 25),
       use_raster = TRUE, raster_quality = 30)
  dev.off()
  
  
  
  # cell_anno_ob用于后续umap，col等
  return(cell_anno_ob)
}

MySubclusterPlot <- function(myinfercnv_obj, cell_anno_ob, output){
  myinfercnv_obj$meta %>% left_join(cell_anno_ob %>% select(cell, subcluster), by = c('barcodes'='cell')) %>%
    filter(str_detect(infercnv_reference, '_observation')) %>%
    mutate(plot_name = str_remove_all(subcluster, '_observation')) %>%
    count(infercnv_group, plot_name) %>%
    ggplot(aes(x = infercnv_group, y = n)) +
    geom_col(aes(fill = plot_name), position = 'fill') +
    scale_fill_manual(values = getPalette_ob(length(unique(cell_anno_ob$subcluster)))) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = 'proportions', x = 'group', title = 'numbers of subcluter cells in group') +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = 'right',
          text = element_text(size = 15),
          legend.text = element_text(size = 6)) -> p
  ggsave(plot = p, filename = str_c(output, '/numbers_of_subcluster_cells_in_group.png'), width = 6, height = 6)
  ggsave(plot = p, filename = str_c(output, '/numbers_of_subcluster_cells_in_group.pdf'), width = 6, height = 6)
  
  MyTsnePlot <- function(cell, cell_anno_ob){
    cell_anno_ob %>%
      filter(str_detect(subcluster, {{cell}})) %>%
      distinct(subcluster) %>% pull(subcluster) %>% length() -> sub_num
    
    if(sub_num > 0){
      print(str_c('plotting ', cell))
      p <- myinfercnv_obj$meta %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(color = 'grey') +
        scale_color_manual(values = getPalette_sample(sub_num)) +
        facet_wrap(~infercnv_group, ncol = 3) +
        labs(title = str_c('subcluters of ', cell)) +
        guides(color = guide_legend(override.aes = list(size = 5))) + 
        geom_point(data = cell_anno_ob %>%
                     filter(str_detect(subcluster, {{cell}}))%>%
                     rename(infercnv_group = group) %>%
                     mutate(plot_name = str_remove_all(subcluster, '_observation')),
                   aes(x = UMAP_1, y = UMAP_2, color = plot_name)) +
        theme_classic() +
        theme(text = element_text(size = 16))
      return(p)
    }
  }
  
  cell_ob <- myinfercnv_obj$meta %>% filter(str_detect(infercnv_reference, '_observation')) %>%
    distinct(cellType) %>% pull(cellType)
  
  for(i in cell_ob){
    p <- MyTsnePlot(cell = i, cell_anno_ob)
    cell <- str_replace_all(i, ' ', '_')
    ggsave(plot = p, filename = str_c(cell, '_subcluster_umap.png'),
           path = output, width = 12, height = 6, limitsize = FALSE)
    ggsave(plot = p, filename = str_c(cell, '_subcluster_umap.pdf'), 
           path = output, width = 12, height = 6, limitsize = FALSE)
  }
  
}


MyCnvScore <- function(myinfercnv_obj, output){
  df <- read.table(myinfercnv_obj$infercnv.observations)
  infer_obj <- readRDS(myinfercnv_obj$infer_obj)
  
  if(length(infer_obj@observation_grouped_cell_indices)>=2){
    meta <- myinfercnv_obj$meta %>% select(barcodes, infercnv_reference)
    df %>% rownames_to_column('genes') %>%
      pivot_longer(cols = -1, names_to = 'barcodes', values_to = 'expr')  %>%
      mutate(barcodes = str_replace_all(barcodes, '\\.', '-')) %>%
      group_by(barcodes) %>%
      summarise(total = sum(expr ^ 2)) %>% 
      left_join(meta, by = 'barcodes') %>%
      ggplot(aes(x = infercnv_reference, y = total)) +
      geom_violin(aes(fill = infercnv_reference)) +
      labs(y = 'CNV scores', x = '') +
      scale_fill_nejm() +
      theme_bw() +
      theme(text = element_text(size = 18), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            legend.position = 'top') -> p
    
  }else{
    print('running k-means')
    kmeans_res <- kmeans(t(as.matrix(df %>% drop_na())), centers = 6)
    results <- as.data.frame(kmeans_res$cluster) %>% 
      rename(cluster = `kmeans_res$cluster`) %>%
      rownames_to_column('barcodes') %>%
      mutate(barcodes = str_replace_all(barcodes, '\\.', '-'))
    
    df %>% rownames_to_column('genes') %>%
      pivot_longer(cols = -1, names_to = 'barcodes', values_to = 'expr')  %>%
      mutate(barcodes = str_replace_all(barcodes, '\\.', '-')) %>%
      group_by(barcodes) %>%
      summarise(total = sum(expr ^ 2)) %>% 
      left_join(results, by = 'barcodes') %>%
      mutate(cluster = str_c('cluster_', as.character(cluster))) %>%
      ggplot(aes(x = cluster, y = total)) +
      geom_violin(aes(fill = cluster)) +
      labs(y = 'CNV scores', x = '') +
      scale_fill_nejm() + 
      theme_bw() +
      theme(text = element_text(size = 18), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            legend.position = 'top') -> p
  }
  
  ggsave(p, filename = 'step3_infercnv_subclusters/score_volin.png')
  ggsave(p, filename = 'step3_infercnv_subclusters/score_volin.pdf')
}


load(opt$myinfercnv_obj)

if(myinfercnv_obj$mode == 'subclusters'){
  MyMkdir('step3_infercnv_subclusters')
  cell_anno_ob <- MyHeatMapPlot(myinfercnv_obj = myinfercnv_obj,
                                output = 'step3_infercnv_subclusters')
  MySubclusterPlot(myinfercnv_obj = myinfercnv_obj,
                   output = 'step3_infercnv_subclusters',
                   cell_anno_ob = cell_anno_ob)

  MyCnvScore(myinfercnv_obj = myinfercnv_obj, output = 'step3_infercnv_subclusters')
}
