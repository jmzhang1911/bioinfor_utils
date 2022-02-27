#!/usr/bin/Rscript

# parameter
suppressMessages(library(getopt))

spec <- matrix(c(
  'qc.rds'     , 'R', 2, 'character', 'Rdata after samples qc, example: single_qc.Rds',
  'od'         , 'o', 2, 'character', 'out dir name',
  'sample'     , 'a', 1, 'character', 'sample name',
  'FDR'        , 'q', 1, 'numeric'  , 'FDR threshold',
  'pvalue'     , 'p', 1, 'numeric'  , 'pvalue', 
  'fold'       , 'f', 1, 'numeric'  , 'fold change threshold',
  'resolution' , 'r', 1, 'numeric'  , 'resolution',
  'min.pct'    , 'm', 1, 'numeric'  , 'gene expression min percent',
  'id'         , 'i', 2, 'character', 'gene id list, example: id_name.list',
  'color'      , 'c', 1, 'character', 'color',
  'size'       , 's', 1, 'numeric'  , 'font size',
  'top.n'      , 'n', 1, 'numeric'  , 'top.n',
  'SCT'        , 'T', 0, 'logical'  , 'use SCT',
  'help'       , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol =5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:AnalysisTsneUmap.R
Version: Version v1.0
Contact: Ranjr <ranjr@biomarker.com.cn> 
Program Date:   2020.9.3 
Description: this program is used to find cluster marker gene ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q()
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$qc.rds) && is.null(opt$od) && is.null(opt$id)) {PrintUsage()}
if(is.null(opt$sample)) {PrintUsage()}
# if(is.null(opt$FDR)) {opt$FDR <- 0.1}
# if(is.null(opt$fold)) {opt$fold <- 2}
if(is.null(opt$resolution)) {opt$resolution <- 0.2}
if(is.null(opt$min.pct)) {opt$min.pct <- 0.25}
if(is.null(opt$color)) {opt$color <- 0.25}
if(is.null(opt$size)) {opt$size <- 10}
if(is.null(opt$top.n)) {opt$top.n <- 10}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

# load packages
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Matrix))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
suppressMessages(library(circlize))
suppressMessages(library(parallel))
suppressMessages(library(ComplexHeatmap))
# library(future)
# plan("multiprocess", workers = 6)
# options(future.globals.maxSize = 1572864000)
# suppressMessages(library(ComplexHeatmap, lib.loc = "/share/nas1/ranjr/packages/3.6"))

time1 <- proc.time()
# parameter set
rds <- opt$qc.rds
od <- opt$od
sample <- opt$sample
FDR <- opt$FDR
pvalue <- opt$pvalue
fold <- opt$fold
res <- opt$resolution
min.pct <- opt$min.pct
ensembleID_symbol <- opt$id
color <- opt$color
size <- opt$size
top.n <- opt$top.n
set.seed(1)

# function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
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

AddId <- function(data, id_list){
  data.re <- inner_join(data, id_list, by = c("gene" = "Symbol")) %>% 
    dplyr::select(EnsembleId, gene, dplyr::everything())
  return(data.re)
}


single.seurat <- readRDS(rds)
ensemble.id.symbol <- read.table(file = ensembleID_symbol, header = FALSE, col.names = c("EnsembleId", "Symbol"))

###zhangjm修改于2021.12.07日，修改逻辑，SCT -> 寻找高可变基因 -> PCA -> 降维聚类
if(!is.null(opt$SCT)){
  single.seurat <- SCTransform(object = single.seurat, do.scale = TRUE, verbose = FALSE, vars.to.regress = c("nCount_RNA"))
  single.seurat <- FindVariableFeatures(single.seurat, selection.method = "vst", nfeatures = 2000)
  my.assay <- 'SCT'
}else{
  single.seurat <- NormalizeData(single.seurat, normalization.method =  "LogNormalize", scale.factor = 10000)
  single.seurat <- FindVariableFeatures(single.seurat, selection.method = "vst", nfeatures = 2000)
  single.seurat <- ScaleData(single.seurat)
  my.assay <- 'RNA'
}

###原始代码：bug：后续PCA结果基于RNA的高可变基因并非SCT的高可变基因
# Findvariabls, Normalized, scaledata and remove batch effect
# single.seurat <- FindVariableFeatures(single.seurat, selection.method = "vst", nfeatures = 2000)
# if(!is.null(opt$SCT)){
  # single.seurat <- SCTransform(object = single.seurat, do.scale = TRUE, verbose = FALSE, vars.to.regress = c("nCount_RNA"))
  # my.assay <- 'SCT'
# }else{
  # single.seurat <- NormalizeData(single.seurat, normalization.method =  "LogNormalize", scale.factor = 10000)
  # single.seurat <- ScaleData(single.seurat)
  # my.assay <- 'RNA'
# }

single.seurat <- RunPCA(single.seurat)
single.seurat <- RunUMAP(single.seurat, dims=1:30)
single.seurat <- RunTSNE(single.seurat, dims=1:30)
single.seurat <- FindNeighbors(single.seurat, dims = 1:30)
single.seurat <- FindClusters(single.seurat, resolution = res) # resolution 
single.seurat[['cellcluster']] <- Idents(single.seurat)[rownames(single.seurat@meta.data)]

# plot
outputdir <- paste(od, "reduction", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}
P.umap <- UMAPPlot(single.seurat, label = TRUE) 
P.tsne <- TSNEPlot(single.seurat, label = TRUE) 
SavePlot(od = outputdir, filename = "umap_all", data = P.umap)
SavePlot(od = outputdir, filename = "tsne_all", data = P.tsne)

pca.projection <- single.seurat@reductions$pca@cell.embeddings
cells <- rownames(pca.projection)
pca.out <- data.frame(Cells = cells, pca.projection)
write.table(pca.out, file = file.path(outputdir, paste(sample, "pca_components.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

tsne.projection <- single.seurat@reductions$tsne@cell.embeddings
cells <- rownames(tsne.projection)
tsne.out <- data.frame(Cells = cells, tsne.projection)
write.table(tsne.out, file = file.path(outputdir, paste(sample, "tsne_components.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

umap.projection <- single.seurat@reductions$umap@cell.embeddings
cells <- rownames(umap.projection)
umap.out <- data.frame(Cells = cells, umap.projection)
write.table(umap.out, file = file.path(outputdir, paste(sample, "umap_components.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# diff marker gene in clusters
logfc.threshold <- log2(fold)
if(!file.exists(file.path(od, "diff_exp_all.Rda"))){
  diff.exp.all <- FindAllMarkers(object = single.seurat, assay = my.assay, slot = "counts", min.pct = 0.25, logfc.threshold = 0)
  
  #2021.12.15
  #diff.exp.all$avg_logFC <- diff.exp.all$avg_logFC/log(2) 
  colnames(diff.exp.all)[2] <- 'avg_logFC'
  
  diff.exp.all <- AddId(diff.exp.all, ensemble.id.symbol)
  save(diff.exp.all, file = file.path(od, "diff_exp_all.Rda"))
}else{
  load(file.path(od, "diff_exp_all.Rda"))
}
print("diff threshold is:")
print(FDR)
print(pvalue)
print(fold)
if(!is.null(FDR)){
  if(!is.null(fold)){
    diff.exp.all.filter <- diff.exp.all %>% filter(p_val_adj < FDR & abs(avg_logFC) > logfc.threshold)
  }else{
    diff.exp.all.filter <- diff.exp.all %>% filter(p_val_adj < FDR)
  }
}
if(!is.null(pvalue)){
  if(!is.null(fold)){
    diff.exp.all.filter <- diff.exp.all %>% filter(p_val < pvalue & abs(avg_logFC) > logfc.threshold)
  }else{
    diff.exp.all.filter <- diff.exp.all %>% filter(p_val < pvalue)
  }
}

# top10 marker gene plot
outputdir <- paste(od, paste("top", top.n, sep = ""), sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}

# plot vln tsane umap
GlobPlotTheme <- function (font.size = size){
  Theme <- theme(title = element_text(size = font.size), 
                   axis.ticks = element_blank(), 
                   axis.line = element_blank(), 
                   axis.title = element_blank(), 
                   axis.text = element_blank())
  return(Theme)
}
#
# parallel::mclapply(sort(unique(Idents(single.seurat))) %>% as.character(), function(x){
#   print(paste(c('doing', x), sep = ' '))
#   top10.genes <- diff.exp.all.filter %>% group_by(cluster) %>% filter(., cluster == x) %>% top_n(., n = top.n, wt = avg_logFC) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),] %>%  .[, "gene", drop = TRUE] %>% head(n = top.n) %>% unique()
#
#   top_gene_num = length(top10.genes)
#   if(top_gene_num < 5){tmp_ncol = top_gene_num}else{tmp_ncol = ceiling(top_gene_num / 2)}
#
#   p.tsne <- FeaturePlot(single.seurat, features = top10.genes, reduction = "tsne", cols = c("lightgrey", "purple"), ncol = tmp_ncol, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
#   outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerTsne", sep = "_")
#   SavePlot(od = outputdir, filename = outfile, data = p.tsne)
#   p.umap <- FeaturePlot(single.seurat, features = top10.genes, reduction = "umap", cols = c("lightgrey", "purple"), ncol = tmp_ncol, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
#   outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerUmap", sep = "_")
#   SavePlot(od = outputdir, filename = outfile, data = p.umap)
#   p.vln <- VlnPlot(single.seurat, features = top10.genes, ncol = 5, pt.size = 0) * theme(title = element_text(size = 7),
#       axis.ticks = element_blank(),
#       axis.title = element_blank(),
#       axis.text.y = element_text(size = 7),
# 	  axis.text.x  = element_text(size = 6, angle = 90))
#   outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerVln", sep = "_")
#   SavePlot(od = outputdir, filename = outfile, data = p.vln)
#   print(paste('---> done', x, sep = ' '))
#   return("done")
# }, mc.cores = length(unique(Idents(single.seurat))))



for(x in c(sort(unique(Idents(single.seurat))) %>% as.character())){
  print(paste('doing', x, sep = ' '))
  top10.genes <- diff.exp.all.filter %>% group_by(cluster) %>% filter(., cluster == x) %>% top_n(., n = top.n, wt = avg_logFC) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),] %>%  .[, "gene", drop = TRUE] %>% head(n = top.n) %>% unique()

  top_gene_num = length(top10.genes)
  if(top_gene_num < 5){tmp_ncol = top_gene_num}else{tmp_ncol = ceiling(top_gene_num / 2)}

  p.tsne <- FeaturePlot(single.seurat, features = top10.genes, reduction = "tsne", cols = c("lightgrey", "purple"), ncol = tmp_ncol, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerTsne", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.tsne)
  p.umap <- FeaturePlot(single.seurat, features = top10.genes, reduction = "umap", cols = c("lightgrey", "purple"), ncol = tmp_ncol, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerUmap", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.umap)
  p.vln <- VlnPlot(single.seurat, features = top10.genes, ncol = 5, pt.size = 0) * theme(title = element_text(size = 7),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 7),
	  axis.text.x  = element_text(size = 6, angle = 90))
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerVln", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.vln)
  print(paste('---> done', x, sep = ' '))
}


# marker gene expression
# marker.gene.avgExp.data <- lapply(sort(unique(Idents(single.seurat))), function(x) {
#   cells <- WhichCells(object = single.seurat, idents = x)
#   t <- as.matrix(apply(single.seurat[[my.assay]]@counts[match(diff.exp.all.filter$gene,rownames(single.seurat[[my.assay]]@counts)), cells], 1, median), ncol = 1)
#   colnames(t) <- x
#   return(t)
# })
avg.data <- AverageExpression(single.seurat)
marker.gene.avgExp <- avg.data$RNA %>% as.data.frame() %>% rownames_to_column("gene")
#marker.gene.avgExp <- do.call(cbind,marker.gene.avgExp.data)
marker.gene.avgExp %<>% as.data.frame %>% AddId(., ensemble.id.symbol)
marker.gene.avgExp <- marker.gene.avgExp[match(unique(diff.exp.all.filter$gene), marker.gene.avgExp$gene),]

# output file
outputdir <- paste(od, "clusterDiff", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}

# plot heatmap
# P.heatmap <- DoHeatmap(single.seurat, features = diff.exp.all.filter %>% 
                       # group_by(cluster) %>% 
                       # top_n(., n = top.n, wt = avg_logFC) %>%
                       # .[, "gene", drop = TRUE] %>% unique(), label = F, draw.lines = F) + theme(axis.text = element_text(size = 3))			
topmarker <- diff.exp.all.filter %>% group_by(cluster) %>% top_n(., n = top.n, wt = avg_logFC) %>% .[, "gene", drop = TRUE] %>% unique()
cluster_info <- sort(single.seurat$seurat_clusters)
exp_data <- GetAssayData(single.seurat, assay = my.assay, slot = "scale.data")
exp_data <- as.matrix(exp_data[intersect(topmarker, rownames(exp_data)), names(cluster_info)])
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = rainbow(length(levels(Idents(single.seurat)))), col = "white"),
                       labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.4, col = "black"))) 
mark_gene <- diff.exp.all.filter %>% group_by(cluster) %>% top_n(., n = 2, wt = avg_logFC) %>% .[, "gene", drop = TRUE] %>% unique()
gene_pos <- match(mark_gene, rownames(exp_data))
grid.draw.Heatmap <- function(x){
  print(x)
}
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = mark_gene, labels_gp = gpar(fontsize = 8)))
P.heatmap <- Heatmap(exp_data,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL,
		col = colorRamp2(c(-5, 1, 5), c("blue","white", "red")),
        heatmap_legend_param = list(
          title = ""
        ))
# plot dotplot					   
P.dotplot <- DotPlot(single.seurat, features = diff.exp.all %>% 
                       group_by(cluster) %>% 
                       top_n(., n = 2, wt = avg_logFC) %>%
                       .[, "gene", drop = TRUE] %>% unique(), dot.min = 0.3) + labs(x = "marker gene", y = "cluster") + theme(axis.text = element_text(size = 10)) + coord_flip()
SavePlot(od = outputdir, filename = "markergene_heatmap", data = P.heatmap)
SavePlot(od = outputdir, filename = "top2_markergene_dotplot", data = P.dotplot)

header <- c("ID", "symbol", "Pvalue", "log2FC", "pct.1", "pct.2", "Qvalue", "Cluster")
# sapply(sort(unique(Idents(single.seurat))) %>% as.character(), function(x){
  # clusterx <- diff.exp.all.filter %>% group_by(cluster) %>% filter(., cluster == x) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),]
  # outfile <- paste(sample, paste("cluster", x, sep = ""),"diff_featuregene.xls", sep = ".")
  # write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # write.table(clusterx, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)	
  # clusterx.all <- diff.exp.all %>% group_by(cluster) %>% filter(., cluster == x) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),]
  # outfile <- paste(paste("cluster", x, sep = ""),"all_featuregene.xls", sep = ".")
  # write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # write.table(clusterx.all, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)
# })

for(x in sort(unique(Idents(single.seurat))) %>% as.character()){
  print(paste('doing', x, sep = ' '))
  clusterx <- diff.exp.all.filter %>% group_by(cluster) %>% filter(., cluster == x) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),]
  outfile <- paste(sample, paste("cluster", x, sep = ""),"diff_featuregene.xls", sep = ".")
  write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(clusterx, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)	
  clusterx.all <- diff.exp.all %>% group_by(cluster) %>% filter(., cluster == x) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),]
  outfile <- paste(paste("cluster", x, sep = ""),"all_featuregene.xls", sep = ".")
  write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(clusterx.all, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)
  print(paste('---> done', x, sep = ' '))
}




outputdir <- od
# write.table(t(header), file = file.path(outputdir, "nofilter_All_clusetr_featuregene.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(diff.exp.all, file = file.path(outputdir, "nofilter_All_clusetr_featuregene.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)
# write.table(t(header), file = file.path(outputdir, "filterd_All_clusetr_markergene.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(diff.exp.all.filter, file = file.path(outputdir, "filterd_All_clusetr_markergene.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)

# statistic marker gene in clusters
up.deg <- diff.exp.all.filter$cluster %>% table() %>% as.matrix() %>% t() %>% as.data.frame()
up.deg.stat <- data.frame(cluster = "deg_number", up.deg)
colnames(up.deg.stat) <- c("cluster", colnames(up.deg))
write.table(up.deg.stat, file = file.path(outputdir, paste(sample, "cluster_deg.stat.xls", sep = ".")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


# marker gene average expression in cluster
header <- c("ID", "symbol", paste("cluster", sort(unique(Idents(single.seurat))), sep = ""))
write.table(t(header), file = file.path(outputdir, paste(sample, "All_cluster_Markergene_avgExp.xls", sep = ".")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(marker.gene.avgExp, file = file.path(outputdir, paste(sample, "All_cluster_Markergene_avgExp.xls", sep = ".")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = T)

cluster.out <- data.frame(rownames(single.seurat@meta.data), single.seurat@meta.data$seurat_clusters)
names(cluster.out) <- c("Barcode", "Cluster")
write.table(cluster.out, file = file.path(outputdir, paste(sample, "clusters.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

counts <- single.seurat@assays$RNA@counts %>% as.data.frame %>% mutate(gene = rownames(.)) %>% AddId(., ensemble.id.symbol)
# scale.data <- single.seurat@assays$RNA@scale.data %>% as.data.frame %>% mutate(gene = rownames(.), .keep = "all") %>% AddId(., ensemble.id.symbol)
write.table(counts, file = file.path(outputdir, paste(sample, "All_cell_counts.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(scale.data, file = file.path(outputdir, "All_cell_expression.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

save(single.seurat, file = file.path(outputdir, "single_analysis.Rda"))
saveRDS(single.seurat, file = file.path(outputdir, "single_seruat.Rds"))

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('Elapsed time: ',run.time[3][[1]],"s"))