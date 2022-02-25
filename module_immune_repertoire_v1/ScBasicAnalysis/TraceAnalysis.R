#!/usr/bin/Rscript

# parameter
suppressMessages(library(getopt))

spec <- matrix(c(
  'seurat'     , 'S', 2, 'character', 'seurat object',
  'assay'      , 'A', 1, 'character', 'setting assays',
  'od'         , 'o', 2, 'character', 'out dir name',
  'sample'     , 'a', 1, 'character', 'sample name',
  'minexp'     , 'e', 1, 'numeric'  , 'min expression',
  'FDR'        , 'f', 1, 'numeric'  , 'FDR change threshold',
  'pvalue'     , 'p', 1, 'numeric'  , 'pvalue threshold',
  'id'         , 'i', 2, 'character', 'gene id list, example: id_name.list',
  'mincell'    , 'C', 1, 'numeric'  , 'min cells',
  'color'      , 'c', 1, 'character', 'color',
  'size'       , 's', 1, 'numeric'  , 'font size',
  'top.n'      , 'n', 1, 'numeric'  , 'top.n',
  'SCT'        , 'T', 0, 'logical'  , 'SCT',
  'anno'       , 'y', 1, 'character', 'annotation column name',
  'cluster'    , 'l', 1, 'character', 'cluster, example: 0,1,2 or B cell',
  'help'       , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol =5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:AnalysisTsneUmap.R
Version: Version v1.0
Contact: Ranjr <ranjr@biomarker.com.cn> 
Program Date:   2020.9.34
Description: this program is used to find cluster marker gene ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q()
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$seurat) & is.null(opt$od)) {PrintUsage()}
if(is.null(opt$sample)) {opt$sample <- "all"}
# if(is.null(opt$FDR)) {opt$FDR <- 0.05}
if(is.null(opt$minexp)) {opt$minexp <- 0.1}
if(is.null(opt$mincell)) {opt$mincell <- 10}
if(is.null(opt$color)) {opt$color <- 0.25}
if(is.null(opt$size)) {opt$size <- 10}
if(is.null(opt$top.n)) {opt$top.n <- 50}
if(!is.null(opt$SCT)) {
  opt$assay <- 'SCT'
}else{
  opt$assay <- 'RNA'
}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

### .libPaths('/share/nas1/ranjr/packages/3.6')
#cmd <- paste("echo ", paste0("'.libPaths('/share/nas1/ranjr/packages/3.6')'"), "> ~/.Rprofile")
#system(cmd, intern = TRUE)

# load packages
suppressMessages(library(Seurat))
suppressMessages(library(proxy))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
### suppressMessages(library(monocle, lib.loc = "/share/nas1/ranjr/packages/3.6"))
suppressMessages(library(monocle))
suppressMessages(library(ggplot2))

# parameter set
seurat.object <- opt$seurat
FDR <- opt$FDR
pvalue <- opt$pvalue
minexp <- opt$minexp
mincell <- opt$mincell
color <- opt$color
size <- opt$size
top.n <- opt$top.n
od <- opt$od
sample <- opt$sample
my.assay <- opt$assay
my.cluster <- opt$cluster
# function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(sample, filename, "png", sep = ".")
  file.pdf <- paste(sample, filename, "pdf", sep = ".")
  tryCatch(
  {ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)},
  error = function(e) {
  tryCatch(
    {
	png(file.path(od, file.png), width = 200, height = 200,res = 300,units = "mm")
    print(data)
    dev.off()
    
    pdf(file.path(od, file.pdf))
    print(data)
    dev.off()
	},
	error = function(e){
	print("我画图又出错了，烦死了！")
	}
  )
  }
)
}

# load(seurat.object)
single.seurat <- readRDS(seurat.object)

# #############################################################原始#######################
# if(!is.null(my.cluster)){
  # single.seurat <- subset(single.seurat, idents = my.cluster)
# }
# # find marker gene
# # single.marker <- FindAllMarkers(single.seurat)
# # single.marker %<>% filter(p_val_adj < FDR)
# # topmarker <- single.marker %>% top_n(., n = top.n, wt = avg_logFC)

# # export seurat to monocle object
# data <- as(as.matrix(single.seurat@assays$RNA@data), "sparseMatrix")
# pd <- new("AnnotatedDataFrame", data = single.seurat@meta.data)
# fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
# fd <- new('AnnotatedDataFrame', data = fData)
# monocle_cds <- newCellDataSet(data,
                              # phenoData = pd,
                              # featureData = fd,
                              # lowerDetectionLimit = minexp,
                              # expressionFamily = negbinomial.size())


# # requred, estimate sizefactor and dispersion
# monocle_cds <- estimateSizeFactors(monocle_cds)
# monocle_cds <- estimateDispersions(monocle_cds)


# # filter low expression genes
# monocle_cds <- detectGenes(monocle_cds, min_expr = minexp)
# expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= mincell))
# monocle_cds <-  monocle_cds[expressed_genes, ]

# disp_table <- dispersionTable(monocle_cds)
# unsup_clustering_genes <- subset(disp_table, mean_expression >= minexp)
# monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)

# # diff_test_res <- differentialGeneTest(monocle_cds,
# #                                      fullModelFormulaStr = "~percent.mt")
# # ordering_genes <- row.names (subset(diff_test_res, qval < 0.05))
# # monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
# varibalegene <- intersect(VariableFeatures(single.seurat), rownames(monocle_cds))
# monocle_cds <- reduceDimension(monocle_cds[varibalegene, ], max_components = 2, 
                               # num_dim = 6, reduction_method = "tSNE")
# monocle_cds <- reduceDimension(monocle_cds[varibalegene, ], max_components = 2,
                               # reduction_method = "DDRTree")
# # monocle_cds <- clusterCells(monocle_cds, num_clusters = length(levels(single.seurat))+ 1)
# monocle_cds <- orderCells(monocle_cds)

# single.marker <- FindAllMarkers(single.seurat, assay = my.assay)
# if(!is.null(opt$FDR)) {single.marker %<>% filter(p_val_adj < FDR)}
# if(!is.null(opt$pvalue)) {single.marker %<>% filter(p_val < pvalue)} 
# print(dim(single.marker))
# # single.marker <- single.marker[intersect(single.marker$gene, rownames(monocle_cds)),] 
# single.marker<- single.marker[single.marker$gene %in% intersect(single.marker$gene, rownames(monocle_cds)),]
# print(head(single.marker))
# topmarker <- single.marker %>% top_n(., n = top.n, wt = avg_logFC) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),] %>% head(n = top.n) 
# print(unique(topmarker$gene))

# select.marker <- single.marker %>% top_n(., n = 5, wt = avg_logFC) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),] %>% .[, "gene", drop = TRUE] %>% head(n = 5) %>% unique()
# sub_cds <- monocle_cds[select.marker, ]

# p.trace.cluster <- plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")
# SavePlot(od = od, filename = "cell_trajectory_cluster", data = p.trace.cluster)
# p.trace.state <- plot_cell_trajectory(monocle_cds, color_by = "State")
# SavePlot(od = od, filename = "cell_trajectory_state", data = p.trace.state)
# p.trace.cluster.marker <- plot_genes_in_pseudotime(sub_cds, color_by = "seurat_clusters")
# SavePlot(od = od, filename = "cell_trajectory_top5_marker_cluster", data = p.trace.cluster.marker)
# p.trace.state.marker <- plot_genes_in_pseudotime(sub_cds, color_by = "State")
# SavePlot(od = od, filename = "cell_trajectory_top5_marker_state", data = p.trace.state.marker)
# # p.trace.heatmap <- plot_pseudotime_heatmap(monocle_cds[unique(topmarker$gene),], show_rownames = TRUE, return_heatmap = TRUE, num_clusters = length(levels(monocle_cds$Cluster)))
# p.trace.heatmap <- plot_pseudotime_heatmap(monocle_cds[unique(topmarker$gene),], show_rownames = TRUE, return_heatmap = TRUE, num_clusters = length(levels(monocle_cds$State)))
# SavePlot(od = od, filename = "cell_trajectory_heatmap", data = p.trace.heatmap)
# save(monocle_cds, file = file.path(od, "cell_trace_monocle.Rda"))

# #cmd <- paste0("sed -i '1,$d' ~/.Rprofile")
# #system(cmd, intern = TRUE)
# #############################################################原始#######################

####修改 ---> zhangjm 2021.11.04 加入判断：若细胞个数大于3w个则跳过此分析过程。
wtf_runner <- function(){

	if(!is.null(my.cluster)){
	  single.seurat <- subset(single.seurat, idents = my.cluster)
	}
	# find marker gene
	# single.marker <- FindAllMarkers(single.seurat)
	# single.marker %<>% filter(p_val_adj < FDR)
	# topmarker <- single.marker %>% top_n(., n = top.n, wt = avg_logFC)

	# export seurat to monocle object
	data <- as(as.matrix(single.seurat@assays$RNA@data), "sparseMatrix")
	pd <- new("AnnotatedDataFrame", data = single.seurat@meta.data)
	fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
	fd <- new('AnnotatedDataFrame', data = fData)
	monocle_cds <- newCellDataSet(data,
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

	# diff_test_res <- differentialGeneTest(monocle_cds,
	#                                      fullModelFormulaStr = "~percent.mt")
	# ordering_genes <- row.names (subset(diff_test_res, qval < 0.05))
	# monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
	varibalegene <- intersect(VariableFeatures(single.seurat), rownames(monocle_cds))
	monocle_cds <- reduceDimension(monocle_cds[varibalegene, ], max_components = 2, 
								   num_dim = 6, reduction_method = "tSNE")
	monocle_cds <- reduceDimension(monocle_cds[varibalegene, ], max_components = 2,
								   reduction_method = "DDRTree")
	# monocle_cds <- clusterCells(monocle_cds, num_clusters = length(levels(single.seurat))+ 1)
	monocle_cds <- orderCells(monocle_cds)

	single.marker <- FindAllMarkers(single.seurat, assay = my.assay)
	if(!is.null(opt$FDR)) {single.marker %<>% filter(p_val_adj < FDR)}
	if(!is.null(opt$pvalue)) {single.marker %<>% filter(p_val < pvalue)} 
	print(dim(single.marker))
	# single.marker <- single.marker[intersect(single.marker$gene, rownames(monocle_cds)),] 
	single.marker<- single.marker[single.marker$gene %in% intersect(single.marker$gene, rownames(monocle_cds)),]
	print(head(single.marker))
	topmarker <- single.marker %>% top_n(., n = top.n, wt = avg_log2FC) %>% .[order(.[,'avg_log2FC'], decreasing = TRUE),] %>% head(n = top.n) 
	print(unique(topmarker$gene))

	select.marker <- single.marker %>% top_n(., n = 5, wt = avg_log2FC) %>% .[order(.[,'avg_log2FC'], decreasing = TRUE),] %>% .[, "gene", drop = TRUE] %>% head(n = 5) %>% unique()
	sub_cds <- monocle_cds[select.marker, ]

	p.trace.cluster <- plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")
	SavePlot(od = od, filename = "cell_trajectory_cluster", data = p.trace.cluster)
	p.trace.state <- plot_cell_trajectory(monocle_cds, color_by = "State")
	SavePlot(od = od, filename = "cell_trajectory_state", data = p.trace.state)
	p.trace.cluster.marker <- plot_genes_in_pseudotime(sub_cds, color_by = "seurat_clusters")
	SavePlot(od = od, filename = "cell_trajectory_top5_marker_cluster", data = p.trace.cluster.marker)
	p.trace.state.marker <- plot_genes_in_pseudotime(sub_cds, color_by = "State")
	SavePlot(od = od, filename = "cell_trajectory_top5_marker_state", data = p.trace.state.marker)
	# p.trace.heatmap <- plot_pseudotime_heatmap(monocle_cds[unique(topmarker$gene),], show_rownames = TRUE, return_heatmap = TRUE, num_clusters = length(levels(monocle_cds$Cluster)))
	p.trace.heatmap <- plot_pseudotime_heatmap(monocle_cds[unique(topmarker$gene),], show_rownames = TRUE, return_heatmap = TRUE, num_clusters = length(levels(monocle_cds$State)))
	SavePlot(od = od, filename = "cell_trajectory_heatmap", data = p.trace.heatmap)
	save(monocle_cds, file = file.path(od, "cell_trace_monocle.Rda"))
}


single.seurat <- readRDS(seurat.object)

if(dim(single.seurat)[2]<= 30000){wtf_runner()}else{print('cells are more than 3w, skipping')}