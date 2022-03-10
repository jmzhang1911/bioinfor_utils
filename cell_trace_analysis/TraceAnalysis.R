#!/usr/bin/env Rscript
# parameter
suppressMessages(library(getopt))

spec <- matrix(c(
  'seurat'     , 'S', 2, 'character', 'seurat object',
  'od'         , 'o', 1, 'character', 'out dir name',
  'cfg'        , 'c', 1, 'character', 'config file',
  'diffdir'   , 'f', 1, 'character', 'deg statistic',
  'minexpminexp'     , 'e', 1, 'numeric'  , 'min expression',
  'mincell'    , 'C', 1, 'numeric'  , 'min cells',
  'top.n'      , 'n', 1, 'numeric'  , 'top.n',
  'perplexity' , 'p', 1, 'numeric'  , 'perplexity of tsne, if the total cells is little, set --perplexity, defalut:30',
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
if(is.null(opt$od)) {opt$od <- "./"}
if(is.null(opt$minexp)) {opt$minexp <- 0.1}
if(is.null(opt$mincell)) {opt$mincell <- 10}
if(is.null(opt$top.n)) {opt$top.n <- 50}
if(is.null(opt$perplexity)) {opt$perplexity <- 30}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

# load packages
suppressMessages(library(Seurat))
suppressMessages(library(proxy))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(monocle))
suppressMessages(library(ggplot2))


# parameter set
para.cfg <- read.delim(file = opt$cfg, sep = "\t", comment.char = "#", row.names = 1, header = FALSE, as.is = TRUE)
para <- para.cfg[,1]
names(para) <- rownames(para.cfg)
sample <- unlist(strsplit(basename(opt$seurat), split = "[.]single_seruat.Rds|[.]integrated_seruat.Rds"))
minexp <- opt$minexp
mincell <- opt$mincell
diffdir <- opt$diffdir
top.n <- opt$top.n
od <- opt$od
dir <- paste(od, paste(sample, "cell_trace", sep = "."), sep = "/")
if(!dir.exists(dir)) {dir.create(dir)}
# function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(sample, filename, "png", sep = ".")
  file.pdf <- paste(sample, filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}


single.seurat <- readRDS(opt$seurat)
# export seurat to monocle object
data <- as(as.matrix(single.seurat@assays$RNA@counts), "sparseMatrix")
pd <- new("AnnotatedDataFrame", data = single.seurat@meta.data)
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

varibalegene <- intersect(VariableFeatures(single.seurat), rownames(monocle_cds))
monocle_cds <- reduceDimension(monocle_cds[varibalegene, ], max_components = 2, perplexity = opt$perplexity,
                               num_dim = 6, reduction_method = "tSNE")
monocle_cds <- reduceDimension(monocle_cds[varibalegene, ], max_components = 2,
                               reduction_method = "DDRTree")
# monocle_cds <- clusterCells(monocle_cds, num_clusters = length(levels(single.seurat))+ 1)
monocle_cds <- orderCells(monocle_cds)
my.data <- pData(monocle_cds)
monocle.data <- my.data %>% rownames_to_column("Cell") %>% select(c("Cell", "sample", "seurat_clusters", "Pseudotime", "State"))
write.table(file = file.path(dir, "Pseudotime_summary.xls"), monocle.data, quote = FALSE, sep = "\t", row.names = FALSE)
p.trace.cluster <- plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")
SavePlot(od = dir, filename = "cell_trajectory_cluster", data = p.trace.cluster)
p.trace.state <- plot_cell_trajectory(monocle_cds, color_by = "State")
SavePlot(od = dir, filename = "cell_trajectory_state", data = p.trace.state)

files <- list.files(diffdir, ".diff_featuregene.xls", full.names = TRUE)
topmarker <- sapply(files, function(file){
  name <- unlist(strsplit(gsub(".diff_featuregene.xls", "", basename(file)), split = "[.]"))[2]
  my.data <- read.table(file = file, sep = "\t", header = T)
  my.data2 <- my.data[match(rownames(monocle_cds), my.data$symbol),]
  select.marker <- my.data2 %>% top_n(., n = 5, wt = log2FC) %>% .[order(.[,'log2FC'], decreasing = TRUE),] %>% .[, "symbol", drop = TRUE] %>% head(n = 5) %>% unique()
  print(select.marker)
if(length(select.marker)>0){ 
  sub_cds <- monocle_cds[select.marker, ]
  print("sub cds")
  p.trace.cluster.marker <- plot_genes_in_pseudotime(sub_cds, color_by = "seurat_clusters")
  SavePlot(od = dir, filename = paste(name, "cell_trajectory_top5_marker_cluster", sep = "."), data = p.trace.cluster.marker)
  p.trace.state.marker <- plot_genes_in_pseudotime(sub_cds, color_by = "State")
  SavePlot(od = dir, filename = paste(name, "cell_trajectory_top5_marker_state", sep = "."), data = p.trace.state.marker)
  return(select.marker)
}
})
topmarker.gene <- unlist(topmarker) %>% as.vector()
print(topmarker.gene)
p.trace.heatmap <- plot_pseudotime_heatmap(monocle_cds[unique(topmarker.gene),], show_rownames = TRUE, return_heatmap = TRUE, num_clusters = length(levels(monocle_cds$State)))
SavePlot(od = dir, filename = "cell_trajectory_heatmap", data = p.trace.heatmap)
#save(monocle_cds, file = file.path(od, "cell_trace_monocle.Rds"))



