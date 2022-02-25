# Name: SingleCellUtils
# Origin data: 2021.08.12
# Updated: 2021.08.13
# Author: zhangjm
# Mail:zhangjm@biomarker.com.cn | jmzhang1911@gmail.com
# Script version: 0.1
# Based on: R4.0.3

library(tidyverse)
library(Seurat)
library(ggplotify)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(RColorBrewer)

MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}

MyPlotsSave <- function(plot, filename='__tmp', output='./', width=10, height=10){
  # Give me a plot, save it as ggplotObj
  MyMkdir(output)
  ggsave(as.ggplot(plot, scale = 0.95), filename = str_c(filename, '.pdf'), path = output, width = width, height = height, limitsize = FALSE)
  width = if_else(width>35, 20, width)
  height = if_else(width>35, 20, height)
  ggsave(as.ggplot(plot, scale = 0.95), filename = str_c(filename, '.png'), path = output, width = width, height = height, dpi = 1200)
}


MyGetGeneExprMatrix <- function(SeuratObj,cell_type_=NA,sample_=NA,gene_symbol_=NA,assay_='RNA',slot_='data'){
  # Give me SeuratObj + celltype + genelist + symbol
  # Return -> gene expression matrix(default assays:RNA slot:counts)
  
  Idents(SeuratObj) <- assay_
  if(is.na(sample_)[1]){sample_ <- SeuratObj@meta.data$sample}else{sample_ <- sample_}
  if(is.na(gene_symbol_)[1]){gene_symbol_ <- rownames(SeuratObj)}else{gene_symbol_ <- gene_symbol_}
  if(is.na(cell_type_)[1]){cell_type_ <- unique(SeuratObj@meta.data$cellType)}else{cell_type_ <- cell_type_}
  
  myfilter <- function(data){dplyr::filter(data, cellType %in% cell_type_, sample %in% sample_, )}
  
  SeuratObj@meta.data %>% rownames_to_column(var = 'barcodes') %>%
    myfilter() %>% pull(barcodes) -> prob
  as.matrix(GetAssayData(SeuratObj, assay = assay_, slot = slot_))[gene_symbol_, prob] %>% as.data.frame() -> df
  return(df)
}


MyHeatMapPlot <- function(in_class_='same',in_matrix_,font_aut_=TRUE,font_size_,
                          main_='__main',legend_title_='__legend_title',
                          heatmap_col_=c('#91CF60','#FFFFBF','#FC8D59'),
                          split_by_,split_col_=pal_nejm(),sample_info_,block_anno_=FALSE){
  # Plotting two kind of heatmaps (default 'same')
  #  type1(same): the x and y are the same, such as gene correlation...
  #  type2(differ): the x and y are different, such as gene expression between different samples
  
  D0_font_aut_ <- function(){
    # find the best fontsize
    # if its not appropriate please add your appropriate parameters below
    # args:
    #   font_aut_: do it or not (default: TRUE), if not font_size_ must be given
    matrix_size <- nrow(in_matrix_)
    if(matrix_size <= 10){
      row_names_gp = 12; column_names_gp = 12
    }else if(matrix_size <= 20){
      row_names_gp = 14; column_names_gp = 14
    }else if(matrix_size <= 30){
      row_names_gp = 6; column_names_gp = 6
    }else if(matrix_size <= 50){
      row_names_gp = 4; column_names_gp = 4
    }else if(matrix_size <= 100){
      row_names_gp = 2; column_names_gp = 2
    }else if(matrix_size <= 600){
      row_names_gp = 1; column_names_gp = 1
    }else{
      row_names_gp = 0.6; column_names_gp = 0.6
    }
    
    font_vector <- c(row_names_gp, column_names_gp)
    names(font_vector) <- c('row_names_gp', 'column_names_gp')
    
    return(font_vector)
  }

  Do_col_aut <- function(){
    # find the default legend scale and color it
    
    p_res <- Heatmap(in_matrix_)
    max_col <- max(p_res@matrix_color_mapping@levels)
    min_col <- min(p_res@matrix_color_mapping@levels)
    col = colorRamp2(breaks = c(min_col, mean(c(max_col, min_col)), max_col),
                     colors = heatmap_col_)
    return(col)
  }
  
  Do_same <- function(){
    # do the type1 heatmap
    # args: 
    #  in_matrix_: input matrix
    #  main_: column title (default: __main)
    #  legend_title_: legend title (default: __legend_title)
    #  heatmap_col_: three colors for heatmap col (default: c('#91CF60','#FFFFBF','#FC8D59'), in order)
    
    Heatmap(in_matrix_, 
            column_names_gp = gpar(fontsize = unit(font_vector[['column_names_gp']], 'npc')),
            row_names_gp = gpar(fontsize = unit(font_vector[['row_names_gp']], 'npc')),
            show_row_dend = F, show_column_dend = F, 
            heatmap_legend_param = list(title = legend_title_,
                                        title_gp = gpar(fontsize = unit(14, 'npc')),
                                        legend_height = unit(0.2, 'npc')),
            column_title = main_, 
            col = Do_col_aut(), 
            column_title_gp = gpar(fontsize = unit(20, 'npc'), 
                                   fontface = 'bold')
            ) -> p_res
    
    return(p_res)
  }
  
  Do_differ <- function(){
    # do the type2 heatmap 
    # args:
    #  in_matrix_: input matrix
    #  sample_info_: for annotation
    #  split_by_: split column by what? must be a column in sample_info_
    #  block_anno_: using annotation or block (default: FALSE), the problem is that the legend overlaped somehow 
    #               when doing column anno in windows 
    #  split_col_: the color of column annotation (default:using pal_nejm()), must be a vector with name
    #              split_col <-c('red','black');names(split_col) <- c('A','B')
    #  main_: column title (default: __main)
    #  legend_title_: legend title (default: __legend_title)
    #  heatmap_col_: three colors for heatmap col (default: c('#91CF60','#FFFFBF','#FC8D59'), in order)
    
    if(is.vector(split_col_)){
      ve = split_col_
    }else{
      ve <- pal_nejm()(length(unique(sample_info_[[split_by_]])))
      names(ve) <- unique(sample_info_[[split_by_]])
    }
    
    if(block_anno_){
      column_annotation <- HeatmapAnnotation(
          cluster = anno_block(gp = gpar(fill = ve), labels = names(ve), 
                               labels_gp = gpar(cex = unit(1.2, 'npc'),
                                                col = 'black',
                                                fontface = 'bold')))
    }else{
      column_annotation <- HeatmapAnnotation(
      cluster = sample_info_[[split_by_]],
      col = list(cluster = ve))
    }

    Heatmap(in_matrix_, 
            column_split = sample_info_[[split_by_]],
            show_row_dend = F, show_column_dend = F, 
            show_column_names = F, show_row_names = T,
            row_names_gp = gpar(fontsize = unit(font_vector[['row_names_gp']], 'npc')),
            top_annotation = column_annotation,
            column_title = main_,
            col = Do_col_aut(),
            heatmap_legend_param = list(title = legend_title_,
                                        title_gp = gpar(fontsize = unit(14, 'npc'),
                                                        fontface = 'bold'),
                                        legend_height = unit(0.3, 'npc'),
                                        title_position = 'leftcenter-rot'),
            column_title_gp = gpar(fontsize = unit(20, 'npc'), 
                                   fontface = 'bold')) -> p_res
    
    return(p_res)
  }
  
  ## main
  if(font_aut_){
    font_vector = D0_font_aut_()
  }else{
    font_vector = font_size_
  }
  
  if(in_class_ == 'same'){
    p <- Do_same()
    return(p)
  }else if(in_class_ == 'differ'){
    p <- Do_differ()
    return(p)
  }
}


MyDimPlot <- function(SeuratObj,reduction){
	getPalette = colorRampPalette(brewer.pal(8, "Set1"))
	DimPlot(SeuratObj, reduction = reduction, label = T) *
		guides(colour = guide_legend(override.aes = list(size = 6))) +
		scale_color_manual(values = getPalette(length(unique(SeuratObj@meta.data$seurat_clusters)))) +
		theme_classic() +
		theme(text = element_text(size = 12, face = 'bold')) -> p
	
	return(p)
}

