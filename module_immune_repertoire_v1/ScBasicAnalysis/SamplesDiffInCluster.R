#!/usr/bin/Rscript

# parameter
### suppressMessages(library(getopt, lib.loc = "/share/nas2/genome/biosoft/R/3.6.1/lib64/R/library/"))
suppressMessages(library(getopt))

###
spec <- matrix(c(
  'rds'        , 'R', 2, 'character', 'Rdata after step1 filter, example: single.Rds',
  'od'         , 'o', 2, 'character', 'out dir name',
  'FDR'        , 'q', 1, 'numeric'  , 'FDR threshold',
  'pvalue'     , 'p', 1, 'numeric'  , 'pvalue', 
  'fold'       , 'f', 1, 'numeric'  , 'fold change threshold',
  'resolution' , 'r', 1, 'numeric'  , 'resolution',
  'min.pct'    , 'm', 1, 'numeric'  , 'gene expression min percent',
  'id'         , 'i', 2, 'character', 'gene id list, example: id_name.list',
  'color'      , 'c', 1, 'character', 'color',
  'size'       , 's', 1, 'numeric'  , 'font size',
  'top.n'      , 'n', 1, 'numeric'  , 'top.n',
  'SCT'        , 'T', 0, 'logical'  , 'SCT',
  'assay'      , 'A', 1, 'character', 'setting assay',
  'scAnno'     , 'a', 1, 'character', 'annotation',
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
if(is.null(opt$rds) && is.null(opt$od) && is.null(opt$id)) {PrintUsage()}
# if(is.null(opt$FDR)) {opt$FDR <- 0.05}
# if(is.null(opt$fold)) {opt$fold <- 2}
if(is.null(opt$resolution)) {opt$resolution <- 0.2}
if(is.null(opt$min.pct)) {opt$min.pct <- 0.25}
if(is.null(opt$color)) {opt$color <- 0.25}
if(is.null(opt$size)) {opt$size <- 10}
if(is.null(opt$top.n)) {opt$top.n <- 10}
#### -> started 2021.10.13 zhangjm 
#### 问题无法递归创建文件夹
# if(!dir.exists(opt$od)) {dir.create(opt$od)} 
#### 修改为：
if(!dir.exists(opt$od)) {dir.create(opt$od, recursive = T)} 


### .libPaths("/share/nas2/genome/biosoft/R/3.6.1/lib64/R/library")
# load packages
suppressMessages(library(SingleCellExperiment))

### suppressMessages(library(Seurat, lib.loc = "/share/nas2/genome/biosoft/R/3.6.1/lib64/R/library"))
suppressMessages(library(Seurat))
###

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Matrix))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
suppressMessages(library(parallel))
suppressMessages(library(circlize))
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 1572864000)

## suppressMessages(library(ComplexHeatmap, lib.loc = "/share/nas1/ranjr/packages/3.6"))
suppressMessages(library(ComplexHeatmap))
### 

set.seed(1)
time1 <- proc.time()
# parameter set
rds <- opt$rds
od <- opt$od
FDR <- opt$FDR
pvalue <- opt$pvalue
fold <- opt$fold
res <- opt$resolution
min.pct <- opt$min.pct
ensembleID_symbol <- opt$id
color <- opt$color
size <- opt$size
top.n <- opt$top.n
if(!is.null(opt$SCT)){
  opt$assay <- 'SCT'
}else{
  opt$assay <- 'RNA'
}
my.assay <- opt$assay
print(paste("use assay is: ", my.assay, sep = ""))
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

grid.draw.Heatmap <- function(x){
  print(x)
}


TopMarkerHeatmapAndDotplot <- function(seuratObject, diff.exp, topN, odir, celltype){
  topmarker <- diff.exp %>% group_by(cluster) %>% top_n(., n = topN, wt = avg_logFC) 
  # cluster_info <- sort(seuratObject$seurat_clusters)
  cluster_info <- sort(Idents(seuratObject))
  # exp_data <- GetAssayData(seuratObject, assay = my.assay, slot = "scale.data")
  exp_data <- GetAssayData(seuratObject, assay = my.assay, slot = "data")
  exp_data <- as.matrix(exp_data[intersect(topmarker$gene, rownames(exp_data)), names(cluster_info)])
  top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = rainbow(length(levels(Idents(seuratObject)))), col = "white"),
                                                     labels = levels(cluster_info), 
                                                     labels_gp = gpar(cex = 0.8, col = "black"))) 
  mark_gene <- topmarker %>% group_by(cluster) %>% top_n(., n = 2, wt = avg_logFC) %>% .[, "gene", drop = TRUE] %>% unique()
  gene_pos <- match(mark_gene, rownames(exp_data))
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
  SavePlot(od = odir, filename = paste("top", topN, "_", "markergene_heatmap", sep = ""), data = P.heatmap)
  # plot dotplot  
  P.dotplot <- DotPlot(seuratObject, assay = my.assay, features = mark_gene, dot.min = 0.3) + labs(x = "marker gene", y = "cluster") + theme(axis.text = element_text(size = 10)) + coord_flip()
  SavePlot(od = odir, filename = "top2_markergene_dotplot", data = P.dotplot)
                       
}  
######################

#load(rda)  
single.integrated <- readRDS(rds)
cell.type <- "seurat_clusters"
if(!is.null(opt$scAnno)){
	cell.type <- opt$scAnno
	Idents(single.integrated) <- single.integrated[[cell.type]]
}
ensemble.id.symbol <- read.table(file = ensembleID_symbol, header = FALSE, col.names = c("EnsembleId", "Symbol"))

# diff marker gene in clusters
###### 修改：添加并行计算
print("diff marker gene in clusters start...")
logfc.threshold <- log2(fold)
if(!file.exists(file.path(od, "diff_exp_all.Rda"))){
  
  ## ---> zhangjm 2021.11.2实现多线程，如果细胞数为6w+原始代码将花费5-6h，修改后20min内可以算完
  # library(BiocParallel)
  # FindMarker.wrapper <- function(x){
    # res = FindMarkers(single.integrated,ident.1= x,  
                      # min.pct=0.25, assay='SCT', slot = "counts", logfc.threshold = 0)
    # res$cluster <- as.character(x)
	# print(paste('---> done the cluster ', as.character(x), sep = ''))
    # return(res)
	# }
	
  # 2021.12.15
  # le <- unique(single.integrated@meta.data$seurat_clusters) %>% as.character() 
  # Markers <- bplapply(le, FindMarker.wrapper, BPPARAM=MulticoreParam(4))
  # diff.exp.all <- do.call(rbind, Markers) %>% rownames_to_column(var = 'gene')
  #diff.exp.all$avg_logFC <- diff.exp.all$avg_logFC/log(2)
  
   # 2021.12.15
  diff.exp.all <- FindAllMarkers(object = single.integrated, assay = my.assay, slot = "counts", min.pct = min.pct, logfc.threshold = 0)
  colnames(diff.exp.all)[2] <- 'avg_logFC'
  diff.exp.all <- AddId(diff.exp.all, ensemble.id.symbol)
  
  save(diff.exp.all, file = file.path(od, "diff_exp_all.Rda"))
  # -------------------------------未探究与FindAllMarker结果是否一致，暂时不用
  

  ###--> 原始方法：
  # diff.exp.all <- FindAllMarkers(object = single.integrated, assay = my.assay, slot = "counts", min.pct = min.pct, logfc.threshold = 0)
  # diff.exp.all$avg_logFC <- diff.exp.all$avg_logFC/log(2)
  # diff.exp.all <- AddId(diff.exp.all, ensemble.id.symbol)
  # save(diff.exp.all, file = file.path(od, "diff_exp_all.Rda"))
  ###--> 原始方法：
  
}else{
  load(file.path(od, "diff_exp_all.Rda"))
}

if(is.null(FDR) && is.null(pvalue)){
  if(!is.null(fold)){
    diff.exp.all.filter <- diff.exp.all %>% filter(abs(avg_logFC) > logfc.threshold)
  }else{
    print("there is no threshold,please check!!!!!")
    q()
  }
}
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
print("diff marker gene in clusters done...")


# top10 marker gene plot
print("top 10 diff marker gene in clusters start...")
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

#####修改zhangjm@2021.11.02，目的跳过没有差异基因的cluster，避免报错
prob <- diff.exp.all.filter %>% count(cluster) %>% filter(n > 0) %>% pull(cluster)
#####原始：sapply(sort(unique(Idents(single.integrated))), function(x){
# sapply(prob, function(x){
  # # diff.features <- diff.exp.all.filter[match(intersect(rownames(single.integrated), diff.exp.all.filter$gene), diff.exp.all.filter$gene),]
  # DefaultAssay(single.integrated) <- my.assay
  # diff.features <- diff.exp.all.filter[diff.exp.all.filter$gene %in% intersect(rownames(single.integrated), diff.exp.all.filter$gene),]
  # deg <- diff.features$cluster %>% table() %>% as.matrix() %>% t() %>% as.data.frame()
  # top10.genes <- diff.features %>% group_by(cluster) %>% filter(., cluster == x) %>% top_n(., n = top.n, wt = avg_logFC) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),] %>% .[, "gene", drop = TRUE] %>% head(n = top.n) %>% unique()
  # len <- length(intersect(top10.genes, rownames(single.integrated)))
  # print(paste("cluster", x,  " intersect gene number: ",len, sep = ""))
  
  # if(len == 0){
    # print(paste('None diff gene in cluster', x, sep=' '))
	# return(paste('None diff gene in cluster', x, sep=' '))
  # }
  # # p.tsne <- FeaturePlot(single.integrated, features = top10.genes, slot = "scale.data", reduction = "tsne", cols = c("lightgrey", "purple"), ncol = 5, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  
  # top_gene_num = length(top10.genes)
  # if(top_gene_num < 5){tmp_ncol = top_gene_num}else{tmp_ncol = ceiling(top_gene_num / 2)}
	
  # p.tsne <- FeaturePlot(single.integrated, features = top10.genes, slot = "data", reduction = "tsne", cols = c("lightgrey", "purple"), ncol = tmp_ncol, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  # outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerTsne", sep = "_")
  # SavePlot(od = outputdir, filename = outfile, data = p.tsne)
  # # p.umap <- FeaturePlot(single.integrated, features = top10.genes, slot = "scale.data", reduction = "umap", cols = c("lightgrey", "purple"), ncol = 5, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  # p.umap <-
# (single.integrated, features = top10.genes, slot = "data", reduction = "umap", cols = c("lightgrey", "purple"), ncol = tmp_ncol, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  # outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerUmap", sep = "_")
  # SavePlot(od = outputdir, filename = outfile, data = p.umap)
  # p.vln <- VlnPlot(single.integrated, features = top10.genes, slot = "data", ncol = 5, pt.size = 0) * theme(title = element_text(size = 7), 
      # axis.ticks = element_blank(), 
      # axis.title = element_blank(), 
      # axis.text.y = element_text(size = 7),
	  # axis.text.x  = element_text(size = 6, angle = 90))
  # outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerVln", sep = "_")
  # SavePlot(od = outputdir, filename = outfile, data = p.vln)
  # finish <- paste("cluster", x, ":", top.n, "gene display is","done")
  # return(finish)
# })

for(x in prob){
  print(paste('doing', x, sep = ' '))
  # diff.features <- diff.exp.all.filter[match(intersect(rownames(single.integrated), diff.exp.all.filter$gene), diff.exp.all.filter$gene),]
  DefaultAssay(single.integrated) <- my.assay
  diff.features <- diff.exp.all.filter[diff.exp.all.filter$gene %in% intersect(rownames(single.integrated), diff.exp.all.filter$gene),]
  deg <- diff.features$cluster %>% table() %>% as.matrix() %>% t() %>% as.data.frame()
  top10.genes <- diff.features %>% group_by(cluster) %>% filter(., cluster == x) %>% top_n(., n = top.n, wt = avg_logFC) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),] %>% .[, "gene", drop = TRUE] %>% head(n = top.n) %>% unique()
  len <- length(intersect(top10.genes, rownames(single.integrated)))
  print(paste("cluster", x,  " intersect gene number: ",len, sep = ""))
  
  if(len == 0){
    print(paste('None diff gene in cluster', x, sep=' '))
	next
  }
  # p.tsne <- FeaturePlot(single.integrated, features = top10.genes, slot = "scale.data", reduction = "tsne", cols = c("lightgrey", "purple"), ncol = 5, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  
  top_gene_num = length(top10.genes)
  if(top_gene_num < 5){tmp_ncol = top_gene_num}else{tmp_ncol = ceiling(top_gene_num / 2)}
	
  p.tsne <- FeaturePlot(single.integrated, features = top10.genes, slot = "data", reduction = "tsne", cols = c("lightgrey", "purple"), ncol = tmp_ncol, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerTsne", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.tsne)
  # p.umap <- FeaturePlot(single.integrated, features = top10.genes, slot = "scale.data", reduction = "umap", cols = c("lightgrey", "purple"), ncol = 5, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  p.umap <- FeaturePlot(single.integrated, features = top10.genes, slot = "data", reduction = "umap", cols = c("lightgrey", "purple"), ncol = tmp_ncol, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerUmap", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.umap)
  p.vln <- VlnPlot(single.integrated, features = top10.genes, slot = "data", ncol = 5, pt.size = 0) * theme(title = element_text(size = 7), 
      axis.ticks = element_blank(), 
      axis.title = element_blank(), 
      axis.text.y = element_text(size = 7),
	  axis.text.x  = element_text(size = 6, angle = 90))
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerVln", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.vln)
  finish <- paste("cluster", x, ":", top.n, "gene display is","done")
  print(finish)
}



print("top 10 diff marker gene in clusters done...")


# output cluster diff analysis file
outputdir <- paste(od, "statistic", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}

# marker gene expression
# 牛逼，lapply里面套apply，apply里面套which！！！这特么是人看的嘛？
marker.gene.avgExp.data <- lapply(sort(unique(Idents(single.integrated))) %>% as.character(), function(x) {
  cells <- WhichCells(object = single.integrated, idents = x)
  #### 修改，可能由于一些定量的id和genelist中的id不匹配导致na
  #### 原始 t <- as.matrix(apply(single.integrated[[my.assay]]@counts[match(diff.exp.all.filter$gene,rownames(single.integrated[[my.assay]]@counts)), cells], 1, median), ncol = 1)
  tmp <- na.omit(match(diff.exp.all.filter$gene,rownames(single.integrated[[my.assay]]@counts)))
  t <- as.matrix(apply(single.integrated[[my.assay]]@counts[tmp, cells], 1, median), ncol = 1)
  colnames(t) <- x
  return(t)
})
marker.gene.avgExp <- do.call(cbind, marker.gene.avgExp.data)
marker.gene.avgExp %<>% as.data.frame %>% mutate(gene = rownames(.)) %>% AddId(., ensemble.id.symbol)
print("marker gene expression in clusters done...")
# marker gene average expression in cluster
header <- c("ID", "symbol", paste("cluster", sort(unique(Idents(single.integrated))), sep = ""))
write.table(t(header), file = file.path(outputdir, "All_cluster_Markergene_avgExp.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(marker.gene.avgExp, file = file.path(outputdir, "All_cluster_Markergene_avgExp.xls"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = T)

# plot heatmap
TopMarkerHeatmapAndDotplot(single.integrated, diff.exp.all.filter, top.n, outputdir, celltype = cell.type)

header <- c("ID", "symbol", "Pvalue", "log2FC", "pct.1", "pct.2", "Qvalue", "Cluster")
for(x in sort(unique(Idents(single.integrated))) %>% as.character() ){
  print(paste('doing', x, sep = ' '))
  clusterx <- diff.exp.all.filter %>% group_by(cluster) %>% filter(., cluster == x) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),]
  outfile <- paste(paste("cluster", x, sep = ""),"diff_featuregene.xls", sep = ".")
  write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(clusterx, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)	
  clusterx.all <- diff.exp.all %>% group_by(cluster) %>% filter(., cluster == x) %>% .[order(.[,'avg_logFC'], decreasing = TRUE),]
  outfile <- paste(paste("cluster", x, sep = ""),"all_featuregene.xls", sep = ".")
  write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(clusterx.all, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)
  print(paste('--> done', x, sep = ' '))
}
# statistic marker gene in clusters
# up.deg.stat <- diff.exp.all.filter$cluster %>% table() %>% as.matrix() %>% t() %>% as.data.frame() %>% mutate(cluster = "up_deg_number", .before = 1)
up.deg <- diff.exp.all.filter$cluster %>% table() %>% as.matrix() %>% t() %>% as.data.frame()
up.deg.stat <- data.frame(cluster = "deg_number", up.deg)
colnames(up.deg.stat) <- c("cluster", colnames(up.deg))
write.table(up.deg.stat, file = file.path(outputdir, "cluster_deg.stat.xls"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('Elapsed time: ',run.time[3][[1]],"s"))
