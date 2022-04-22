library(optparse)

option_list <- list(
  make_option(c('-s', '--seurat_obj'), type = 'character', help = 'seurat object'),
  make_option(c('-c', '--cell_type'), type = 'character', help = 'colname of celltype', default = 'cellType'),
  make_option(c('-o', '--output'), type = 'character', help = 'output', default = 'cell_trace_results')
)

opt <- parse_args(OptionParser(option_list = option_list))
library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
library(Seurat)
library(monocle)


MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}

monocle_trace <- function(seob, dir, cell_type){
  print('doing monocle cell trajectory')
  
  if(file.exists(file.path(dir, 'monocle_cds.RData'))){
    print('loading monocle_cds.RData')
    load(file.path(dir, 'monocle_cds.RData'))
    return(monocle_cds)
  }
  
  MyMkdir(dir)
  minexp <- 0.1
  mincell <- 10
  top.n <- 50
  perplexity <- 30
  
  single.seurat <- readRDS(seob)
  single.seurat@meta.data %>% filter(!is.na(cell_type)) %>% rownames() -> bar
  single.seurat <- single.seurat[, bar]
  
  # 生成用于报告的配置文件
  config <- data.frame(cell_number = dim(single.seurat)[2],
                       cell_type_number = table(single.seurat@meta.data[cell_type]) %>% length())
  write.table(t(config), file = 'config.txt', row.names = T, quote = F, sep = '\t', col.names = F)

  
  if(!cell_type %in% colnames(single.seurat@meta.data)){
    stop('wrong celltype colname')
  }
  
  VariableFeatures <- VariableFeatures(single.seurat)
  # export seurat to monocle object
  data <- as(as.matrix(single.seurat@assays$RNA@counts), "sparseMatrix")
  pd <- new("AnnotatedDataFrame", data = single.seurat@meta.data)
  
  
  
  rm(single.seurat);gc()
  
  fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  monocle_cds <- newCellDataSet(as(data, "sparseMatrix"),
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = minexp,
                                expressionFamily = negbinomial.size())
  
  # requred, estimate sizefactor and dispersion
  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)
  
  # filter low expression genes
  monocle_cds <- detectGenes(monocle_cds, min_expr = minexp)
  expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= mincell))
  monocle_cds <-  monocle_cds[expressed_genes, ]
  
  disp_table <- dispersionTable(monocle_cds)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= minexp)
  monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
  
  varibalegene <- intersect(VariableFeatures, rownames(monocle_cds))
  monocle_cds <- reduceDimension(monocle_cds[varibalegene, ], max_components = 2, perplexity = perplexity,
                                 num_dim = 6, reduction_method = 'tSNE')
  monocle_cds <- reduceDimension(monocle_cds[varibalegene, ], max_components = 2,
                                 reduction_method = "DDRTree")
  # monocle_cds <- clusterCells(monocle_cds, num_clusters = length(levels(single.seurat))+ 1)
  monocle_cds <- orderCells(monocle_cds)
  my.data <- pData(monocle_cds)
  monocle.data <- my.data %>% 
    rownames_to_column("Cell") %>% 
    select(c("Cell", "sample", all_of(cell_type), "Pseudotime", "State"))

  p.trace.cluster <- plot_cell_trajectory(monocle_cds, color_by = cell_type)
  SavePlot(od = dir, filename = "cell_trajectory_cluster", data = p.trace.cluster)
  p.trace.state <- plot_cell_trajectory(monocle_cds, color_by = "State")
  SavePlot(od = dir, filename = "cell_trajectory_state", data = p.trace.state)
  
  
  BEAM_res <- BEAM(monocle_cds, branch_point = 1, cores = 4)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  
  pdf(file = file.path(dir, 'plot_genes_branched_heatmap.pdf'))
  plot_genes_branched_heatmap(monocle_cds[row.names(subset(BEAM_res, qval < 1e-4)),],
                              cores = 4,
                              use_gene_short_name = T,
                              show_rownames = F)
  dev.off()
  
  png(filename = file.path(dir, 'plot_genes_branched_heatmap.png'), res = 200, width = 1200, height = 1200)
  plot_genes_branched_heatmap(monocle_cds[row.names(subset(BEAM_res, qval < 1e-4)),],
                              cores = 4,
                              use_gene_short_name = T,
                              show_rownames = F)
  dev.off()
  
  print('done monocle cell trajectory')
  save(monocle_cds, varibalegene, file = file.path(dir, 'monocle_cds.RData'))
  write.table(file = file.path(dir, "Pseudotime_summary.xls"), monocle.data, quote = FALSE, sep = "\t", row.names = FALSE)
  
  return(monocle_cds)
}


plot_marker <- function(seob, monocle_cds, dir, cell_type){
  print('doing plotting')
  seob <- readRDS(seob)
  seob <-  SetIdent(seob, value = cell_type)
  
  seob@meta.data %>% filter(!is.na(cell_type)) %>% rownames() -> bar
  seob <- seob[, bar]
  
  if(!file.exists(file.path(dir, 'diff.exp.all_origin.RData'))){
    print('doing FindAllMarkers')
    library(future)
    
    plan("multiprocess", workers = 10)
    options(future.globals.maxSize = 1572864000) 
    diff.exp.all_origin <- FindAllMarkers(object = seob, 
                                          min.pct = 0.25, 
                                          logfc.threshold = 0)
    save(diff.exp.all_origin, file = file.path(dir, 'diff.exp.all_origin.RData'))
  }else{
    print('loading FindAllMarkers')
    load(file.path(dir, 'diff.exp.all_origin.RData'))
  }
  
  diff.exp.all_origin %>% as_tibble() %>%
    filter(p_val_adj <= 0.1) %>%
    group_by(cluster) %>% slice_max(avg_log2FC, n = 5) -> markers
  
  for(cluster_tmp in unique(markers$cluster)){
    name <- str_replace_all(cluster_tmp, ';| |,|-', '_')
    markers %>% filter(cluster == cluster_tmp) %>%
      pull(gene) %>% unique() -> select.marker
    
    sub_cds <- monocle_cds[select.marker, ]
    p.trace.cluster.marker <- plot_genes_in_pseudotime(sub_cds, color_by = cell_type)
    SavePlot(od = dir, filename = paste(name, "cell_trajectory_top5_marker_cluster", sep = "."), data = p.trace.cluster.marker)
    p.trace.state.marker <- plot_genes_in_pseudotime(sub_cds, color_by = "State")
    SavePlot(od = dir, filename = paste(name, "cell_trajectory_top5_marker_state", sep = "."), data = p.trace.state.marker)
  }
  
  p.trace.heatmap <- plot_pseudotime_heatmap(
    monocle_cds[unique(unique(markers$gene)),], 
    show_rownames = TRUE, 
    return_heatmap = TRUE, 
    num_clusters = length(levels(monocle_cds$State)))
  SavePlot(od = dir, filename = "cell_trajectory_heatmap", data = p.trace.heatmap)
  
  save(markers, file = file.path(dir, 'markers.RData'))
}

monocle_cds <- monocle_trace(seob = opt$seurat_obj, dir = opt$output, cell_type = opt$cell_type)
plot_marker(seob = opt$seurat_obj,
            monocle_cds = monocle_cds,
            dir = opt$output,
            cell_type = opt$cell_type)
