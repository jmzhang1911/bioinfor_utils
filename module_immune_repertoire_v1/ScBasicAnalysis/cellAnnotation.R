#!/usr/bin/Rscript
suppressMessages(library(getopt))

spec <- matrix(c(
  'object', 'r', 2, 'character', 'seurat Rds oject',
  'ref'   , 'f', 2, 'character', 'reference database',
  'od'    , 'o', 1, 'character', 'outdir name',
  'sample', 's', 2, 'character', 'sample name',
  'help'  , 'h', 0, 'logical'  , 'print help information'
), byrow = T, ncol = 5)

opt <- getopt(spec)
PrintUsage <- function(){
  cat("ProgramName:cellAnnotation.R
Version: Version v1.0
Contact: Ranjr <ranjr@biomarker.com.cn> 
Program Date:   2020.10.16
Description: this program is used to cell annotation ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q()
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$object) & is.null(opt$ref)) {PrintUsage()}
if(is.null(opt$od)) {opt$od <- "./"}
if(!dir.exists(opt$od)) {dir.create(opt$od)}
### print(.libPaths())
suppressMessages(library(Seurat))
suppressMessages(library(SingleR))
suppressMessages(library(ggplot2))
suppressMessages(library(scater))

### suppressMessages(library(scran,lib.loc="/share/nas1/ranjr/packages/3.6"))
suppressMessages(library(scran))
###

seurat.object <- opt$object
ref.data <- opt$ref
od <- opt$od
samples <- opt$sample
seurat.data <- readRDS(seurat.object)
ref.data <- readRDS(ref.data)
print("read data")
test.data <- GetAssayData(object = seurat.data, slot = "data")
com <- intersect(rownames(test.data), rownames(ref.data))

by.t <- scran::pairwiseTTests(assay(ref.data, 1), ref.data$label.main)
markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=100)
# pred.data <- SingleR(test = test.data, ref = ref.data, labels = ref.data$label.main, genes = markers)
pred.data <- SingleR(test = test.data, ref = ref.data, labels = ref.data$label.main,  method =  "cluster", clusters = seurat.data@active.ident)

score <- pred.data$scores
rownames(score) <- paste("cluster", rownames(pred.data), sep = "")
header <- c("clusterID", colnames(score))
write.table(t(header), file = file.path(od, paste0(samples,".cluster_annotation_scores.xls")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(score, file = file.path(od, paste0(samples,".cluster_annotation_scores.xls")), sep = "\t", row.names = TRUE, quote = FALSE, col.names = FALSE, append = TRUE)

result <- data.frame(clusterID = rownames(pred.data), cellAnnotationType = pred.data$pruned.labels)
write.table(result, file = file.path(od,paste0(samples,".cluster_annotation_result.xls")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# p.out.dis <- plotScoreDistribution(results = pred.data)
# ggsave(file = file.path(od, "cell_annotation_ScoreDistribution.png"), p.out.dis)
# p.out <- plotScoreHeatmap(pred.data)
p.out <- plotScoreHeatmap(pred.data, fontsize.row = 9,show_colnames = T)
ggsave(file = file.path(od, paste0(samples,".cluster_annotation_heatmap.png")), p.out)
ggsave(file = file.path(od, paste0(samples,".cluster_annotation_heatmap.pdf")), p.out)

seurat.data[['cellType']] <- "NA"
for(i in 1:nrow(result)){
  seurat.data$cellType[which(seurat.data$seurat_clusters == result[i, 'clusterID'])] <- as.character(result[i, 'cellAnnotationType'])
}

print(unique(seurat.data$cellType))

save(seurat.data, file = file.path(od, "singleCell_annotation.Rda"))
saveRDS(seurat.data, file = file.path(od, "singleCell_annotation.Rds"))

p.out <- DimPlot(seurat.data, group.by="cellType", label=T, reduction='tsne', label.size = 3)
ggsave(file = file.path(od, paste0(samples,".cluster_annotation_tsne.png")), p.out, width = 6, height = 4, scale = 1.3)
ggsave(file = file.path(od, paste0(samples,".cluster_annotation_tsne.pdf")), p.out, width = 6, height = 4, scale = 1.3)
p.out <- DimPlot(seurat.data, group.by="cellType", label=T, reduction='umap', label.size = 3)
ggsave(file = file.path(od, paste0(samples,".cluster_annotation_umap.png")), p.out, width = 6, height = 4, scale = 1.3)
ggsave(file = file.path(od, paste0(samples,".cluster_annotation_umap.pdf")), p.out, width = 6, height = 4, scale = 1.3)
























