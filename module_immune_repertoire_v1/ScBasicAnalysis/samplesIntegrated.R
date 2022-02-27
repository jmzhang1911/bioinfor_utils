#!/usr/bin/Rscript

# parameter
suppressMessages(library(getopt))

spec <- matrix(c(
  'rda'        , 'd', 1, 'character', 'Rda after step1 filter, for example: upload.Rdata',
  'rds'        , 'R', 1, 'character', 'Rds after step1 filter, for example: /path/sample0/human_sce_qc.Rds,/path/sample1/human_sce_qc.Rds, with --mix',
  'od'         , 'o', 1, 'character', 'out dir name',
  'resolution' , 'r', 1, 'numeric'  , 'resolution',
  'mix'        , 'm', 0, 'logical'  , '--mix|-m, mixed data',
  'id'         , 'i', 2, 'character', 'gene id list, example: id_name.list',
  'SCT'        , 'T', 0, 'logical'  , 'SCT',
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
if(is.null(opt$rda) && is.null(opt$rds)) {PrintUsage()}
if(is.null(opt$od) && is.null(opt$id)) {PrintUsage()}
if(is.null(opt$FDR)) {opt$FDR <- 0.05}
if(is.null(opt$fold)) {opt$fold <- 2}
if(is.null(opt$resolution)) {opt$resolution <- 0.1}
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
suppressMessages(library(parallel))
suppressMessages(library(circlize))
suppressMessages(library(RColorBrewer))
# library(future)
# plan("multiprocess", workers = 6)
# options(future.globals.maxSize = 1572864000)

time1 <- proc.time()
# parameter set
rda <- opt$rda
rds <- opt$rds
od <- opt$od
res <- opt$resolution
ensembleID_symbol <- opt$id
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

if(!is.null(opt$mix)){ 
  if(is.null(opt$rds)){
    print("Error:when you set --mix, --rds is must, please check input!!!!")
    q()
  }
  if(!is.null(opt$rda)){
    print("Error:when you set --mix, --rda must be null, please check input!!!!")
    q()
  }
  RDS <- strsplit(rds, split = ",") %>% unlist()
  singleList <- sapply(RDS, function(x){
    object <- readRDS(x)
    sample <- str_split(basename(dirname(x)), '\\.')[[1]][1]
    object[['sample']] <- sample
    return(object)
  })
  names(singleList) <- basename(dirname(RDS))
}

if(!is.null(opt$rda) & is.null(opt$mix)){
  load(rda)
}

sample <- names(singleList)
ensemble.id.symbol <- read.table(file = ensembleID_symbol, header = FALSE, col.names = c("EnsembleId", "Symbol"))
# Findvariabls, Normalized, scaledata and remove batch effect
if(!is.null(opt$SCT)){
  for (i in 1:length(singleList)){
    singleList[[i]] <- SCTransform(object = singleList[[i]], do.scale = TRUE, verbose = FALSE, vars.to.regress = c("nCount_RNA"))
  }
  pancreas.features <- SelectIntegrationFeatures(object.list = singleList, nfeatures = 3000)
  pancreas.list <- PrepSCTIntegration(object.list = singleList, anchor.features = pancreas.features)
  print("Done")
  single.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = pancreas.features, dims = 1:30, normalization.method = "SCT")
  single.integrated <- IntegrateData(anchorset = single.anchors, dims = 1:30)
  DefaultAssay(single.integrated) = 'integrated' # zhangjm
  print("data integrated done...")
}else{
  merge.y <- function(object.list) {
    y = c(object.list[[2]])
    if (length(object.list) > 2L) {
      for (i in 3L:length(object.list)) {
        y <- c(y, object.list[[i]])
      }
    }
    return(y)
  }
  single.integrated <- merge(x = singleList[[1]], y = merge.y(singleList), add.cell.ids = sample)
  single.integrated <- FindVariableFeatures(single.integrated, selection.method = "vst", nfeatures = 2000)
  DefaultAssay(single.integrated) = 'integrated' # zhangjm
}

# VlnPlot scatter 
outputdir <- paste(od, "base", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}
single.integrated@meta.data$nUMI <- single.integrated@meta.data$nCount_RNA 
single.integrated@meta.data$nGene <- single.integrated@meta.data$nFeature_RNA 
P.vlnplot <- VlnPlot(single.integrated,
                     features = c("nGene", "nUMI", "percent.mt"), 
                     ncol = 3, pt.size = 0.1, group.by = "sample") * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10))
SavePlot(od = outputdir, filename = paste("All", "vln", sep = "_"), data = P.vlnplot)

P1.scatter <- FeatureScatter(single.integrated, feature1 = "nUMI", feature2 = "percent.mt", pt.size = 0.5, group.by = "sample") + theme(legend.position = "none", axis.text = element_text(size = 10))
P2.scatter <- FeatureScatter(single.integrated, feature1 = "nUMI", feature2 = "nGene", pt.size = 0.5, group.by = "sample") + theme(axis.text = element_text(size = 10))
P.scatter <- P1.scatter + P2.scatter
SavePlot(od = outputdir, filename = paste("All", "scatter", sep = "_"), data = P.scatter)
print("VlnPlot and scatter done...")

single.integrated <- ScaleData(single.integrated, verbose = FALSE)
single.integrated <- RunPCA(single.integrated, features = VariableFeatures(single.integrated))

# PCA plot
outputdir <- paste(od, "reduction", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}
single.integrated <- JackStraw(single.integrated, num.replicate = 100)
single.integrated <- ScoreJackStraw(single.integrated, dims = 1:20)
# P.Jack <- JackStrawPlot(single.integrated, dims = 1:20)
P.Elbow <- ElbowPlot(single.integrated, ndims = 30)
P.pca <- DimPlot(single.integrated, reduction  = "pca", group.by = "sample")
# SavePlot(od = outputdir, filename = "JackStrawPlot", data = P.Jack)
SavePlot(od = outputdir, filename = "ElbowPlot", data = P.Elbow)
SavePlot(od = outputdir, filename = "PCA_dim", data = P.pca)

###### Started ---> zhangjm修改2021.10.12
### 问题：使用JackStraw分析之后进行UMAP分析会导致UAMP使用的assay为RNA
### 修改方案：指定UAMP的assay以及dim，不再使用features，dim的值其实影响不大，由于SCT对PCA较敏感因此使用1:30
### 源码：
# single.integrated <- RunUMAP(single.integrated, features = VariableFeatures(single.integrated))
# single.integrated <- RunTSNE(single.integrated, features = VariableFeatures(single.integrated))
single.integrated <- RunUMAP(single.integrated, assay = 'integrated', dims = 1:30)
single.integrated <- RunTSNE(single.integrated, assay = 'integrated', dims = 1:30)
#
single.integrated <- FindNeighbors(single.integrated, dims = 1:30)
single.integrated <- FindClusters(single.integrated, resolution = res) # resolution 
single.integrated[['cellcluster']] <- Idents(single.integrated)[rownames(single.integrated@meta.data)]
###
###### End ---> zhangjm修改2021.10.12


# UMAP/TSNE plot
P.umap.all <- UMAPPlot(single.integrated, label = TRUE) 
P.umap.split <- UMAPPlot(single.integrated, split.by = "sample", label = TRUE)
P.umap.sample <- UMAPPlot(single.integrated, group.by = "sample", label = TRUE)
P.tsne.all <- TSNEPlot(single.integrated, label = TRUE) 
P.tsne.split <- TSNEPlot(single.integrated, split.by = "sample", label = TRUE) 
P.tsne.sample <- TSNEPlot(single.integrated, group.by = "sample") 
SavePlot(od = outputdir, filename = "umap_all", data = P.umap.all)
SavePlot(od = outputdir, filename = "umap_splt", data = P.umap.split)
SavePlot(od = outputdir, filename = "umap_sample", data = P.umap.sample)
SavePlot(od = outputdir, filename = "tsne_all", data = P.tsne.all)
SavePlot(od = outputdir, filename = "tsne_splt", data = P.tsne.split)
SavePlot(od = outputdir, filename = "tsne_sample", data = P.tsne.sample)

pca.projection <- single.integrated@reductions$pca@cell.embeddings
cells <- rownames(pca.projection)
pca.out <- data.frame(Cells = cells, pca.projection)
write.table(pca.out, file = file.path(outputdir, "pca_components.xls"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

tsne.projection <- single.integrated@reductions$tsne@cell.embeddings
cells <- rownames(tsne.projection)
tsne.out <- data.frame(Cells = cells, tsne.projection)
write.table(tsne.out, file = file.path(outputdir, "tsne_components.xls"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

umap.projection <- single.integrated@reductions$umap@cell.embeddings
cells <- rownames(umap.projection)
umap.out <- data.frame(Cells = cells, umap.projection)
write.table(umap.out, file = file.path(outputdir, "umap_components.xls"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
print("reductions done...")

# output other file
outputdir <- od
cluster.out <- data.frame(rownames(single.integrated@meta.data), single.integrated@meta.data$seurat_clusters)
names(cluster.out) <- c("Barcode", "Cluster")
write.table(cluster.out, file = file.path(outputdir, paste("All", "ncells_clusters.xls", sep = "_")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

avg.data <- AverageExpression(single.integrated)
marker.gene.avgExp <- avg.data$RNA %>% as.data.frame() %>% rownames_to_column("gene")
print(dim(marker.gene.avgExp))
marker.gene.avgExp <- marker.gene.avgExp %>% as.data.frame %>% AddId(., ensemble.id.symbol)
print(head(marker.gene.avgExp))
header <- c("ID", "symbol", paste("cluster", sort(unique(Idents(single.integrated))), sep = ""))
write.table(t(header), file = file.path(outputdir, "All_cluster_Markergene_avgExp.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(marker.gene.avgExp, file = file.path(outputdir, "All_cluster_Markergene_avgExp.xls"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = T)




# clusters
# res <- 0.2
if(!is.null(opt$SCT)) {
  col.nam <- paste("integrated_snn_res", res, sep = ".")
}else{
  col.nam <- paste("RNA_snn_res", res, sep = ".")
}
stat.table <- table(single.integrated@meta.data[, col.nam], single.integrated@meta.data$sample)
sample.cluster <- matrix(stat.table, ncol = length(unique(single.integrated@meta.data$sample)))
sample.cluster.out <- data.frame(rownames(stat.table), sample.cluster)
names(sample.cluster.out) <- c("cluster", colnames(stat.table))
write.table(sample.cluster.out, file = file.path(outputdir, paste("All", "ncells_clusters.stat.xls", sep = "_")), sep = "\t", row.names = F, quote = FALSE)
if(length(unique(single.integrated@meta.data$sample)) <= 8){
  P.bar <- ggplot(data = single.integrated@meta.data, aes(x = seurat_clusters, fill = sample)) + geom_bar(position = "dodge") + 
    scale_fill_manual(values = brewer.pal(8, "Set2")[1:length(unique(single.integrated@meta.data$sample))]) + 
    theme_bw() + theme(panel.grid=element_blank(),panel.border =element_blank(),axis.line=element_line(size=0.5,colour="black"))
}else{
  P.bar <- ggplot(data = single.integrated@meta.data, aes(x = seurat_clusters, fill = sample)) + geom_bar(position = "dodge") + 
  scale_fill_manual(values = rainbow(length(unique(single.integrated@meta.data$sample)), alpha = 0.8)) +
  theme_bw() + theme(panel.grid=element_blank(),panel.border =element_blank(),axis.line=element_line(size=0.5,colour="black"))
}
SavePlot(od = outputdir, filename = "sample_barplot", data = P.bar)

save(single.integrated, file = file.path(outputdir, "single_analysis.Rda"))
saveRDS(single.integrated, file = file.path(outputdir, "single_seruat.Rds"))

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('Elapsed time: ',run.time[3][[1]],"s"))
