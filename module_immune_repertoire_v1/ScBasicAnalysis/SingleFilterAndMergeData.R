#!/usr/bin/Rscript

# parameter
suppressMessages(library(getopt))
spec <- matrix(c(
  'input.dirs'          , 'i', 2, 'character', 'input dir, Comma separated',
  'od'                  , 'o', 2, 'character', 'output dir name',
  'nUMI.min'            , 'u', 1, 'integer'  , 'min UMI number, default 100',
  'nGene.min'           , 'g', 1, 'integer'  , 'min gene number, default 500',
  'nGene.max'           , 'm', 1, 'integer'  , 'max gene number, default none',
  'percent.mt.max'      , 'P', 1, 'numeric'  , 'mt percent, default 0.2',
  'min.cells'           , 'c', 1, 'integer'  , 'gene expression\'s min cells',
  'gene'                , 'f', 1, 'numeric'  , 'the column of gene,default:2',
  'help'                , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:filterAndMergeData.R
Version: Version v1.0
Contact: Ranjr <ranjr@biomarker.com.cn> 
Program Date:   2020.9.2 
Description: this program is used to filter single cell data ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$input.dirs) & is.null(opt$od)) {PrintUsage()}
if(is.null(opt$nUMI.min)) {opt$nUMI.min <- 100}
if(is.null(opt$nGene.min)) {opt$nGene.min <- 500}
if(is.null(opt$percent.mt.max)) {opt$percent.mt.max <- 0.2}
if(is.null(opt$min.cells)) {opt$min.cells <- 10}
if(is.null(opt$gene)) {opt$gene <- 2}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

### .libPaths("/share/nas2/genome/biosoft/R/3.6.1/lib64/R/library")
# load packages
suppressMessages(library(SingleCellExperiment))
### suppressMessages(library(Seurat, lib.loc="/share/nas2/genome/biosoft/R/3.6.1/lib64/R/library"))
suppressMessages(library(Seurat))

suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
suppressMessages(library(RCurl))
suppressMessages(library(parallel))

# parameter set
nUMI.min <- opt$nUMI.min
nGene.min <- opt$nGene.min
nGene.max <- opt$nGene.max
pct.mt.max <- opt$percent.mt.max
cells.min <- opt$min.cells
genecol <- opt$gene
od <- opt$od
object <- strsplit(opt$input.dirs, split = ",")[[1]]
my.vec <- unlist(strsplit(object, split = "/outs"))
#sample <- basename(my.vec[seq(1, length(my.vec), 2)])
sample <- str_split(basename(my.vec),'\\.')[[1]][1]

print(sample)
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

# parallel computing
outputdir <- paste(od, "singleSample", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}
print(length(object))

g.file <- paste(object[1], "filtered_feature_bc_matrix/features.tsv.gz", sep = "/")
sym <- read.csv(gzfile(g.file), sep = "\t", header = F)
colnames(sym) <- c("ID", "Symbol", "other")
rownames(sym) <- sym$ID
sym <- sym[,1:2]
write.table(sym, file = file.path(od, "symbol.list"), quote = F, row.names = F, col.names = F,  sep = "\t")

singleList <- parallel::mclapply(object, function(x) {
# singleList <- sapply(object, function(x) {
  filter.dir <- paste(x, "filtered_feature_bc_matrix", sep = "/")
  single.matrix <- Read10X(data.dir = filter.dir, gene.column = genecol)
  sample.name <- basename(strsplit(x, split = "/outs")[[1]][1])
  singledir <- paste(outputdir, str_split(sample.name,'\\.')[[1]][1], sep = "/")
  if(!dir.exists(singledir)) {dir.create(singledir)}
  seurat.object <- CreateSeuratObject(counts = single.matrix, project = sample.name, min.cells = cells.min)
  seurat.object[['sample']] <- str_split(sample.name,'\\.')[[1]][1]
  seurat.object$"percent.mt" <- PercentageFeatureSet(object = seurat.object, pattern = "^MT-|^mt-")
  seurat.object$"percent.mt" <- seurat.object$"percent.mt" /100
  
  seurat.object@meta.data$nUMI <- seurat.object@meta.data$nCount_RNA 
  seurat.object@meta.data$nGene <- seurat.object@meta.data$nFeature_RNA 
  P.vlnplot.pre <- VlnPlot(seurat.object,
                    features = c("nGene", "nUMI", "percent.mt"), 
                    ncol = 3, pt.size = 0.1) * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10))
  SavePlot(od = singledir, filename = paste(sample.name, "nofilter.vln", sep = "_"), data = P.vlnplot.pre)
  
  P1.scatter.pre <- FeatureScatter(seurat.object, feature1 = "nUMI", feature2 = "percent.mt", pt.size = 0.5) + theme(legend.position = "none", axis.text = element_text(size = 10))
  P2.scatter.pre <- FeatureScatter(seurat.object, feature1 = "nUMI", feature2 = "nGene", pt.size = 0.5) + theme(axis.text = element_text(size = 10))
  P.scatter.pre <- P1.scatter.pre + P2.scatter.pre
  SavePlot(od = singledir, filename = paste(sample.name, "nofilter.scatter", sep = "_"), data = P.scatter.pre)
  if(!is.null(opt$nGene.max)){
    single.seurat <- subset(x = seurat.object, subset= ((nFeature_RNA >= nGene.min) & (nFeature_RNA <= nGene.max) & (percent.mt < pct.mt.max) & (nCount_RNA >= nUMI.min)))
  }else{
    single.seurat <- subset(x = seurat.object, subset= ((nFeature_RNA >= nGene.min) & (percent.mt < pct.mt.max) & (nCount_RNA >= nUMI.min)))
  }
  P.vlnplot <- VlnPlot(single.seurat,
                    features = c("nGene", "nUMI", "percent.mt"), 
                    ncol = 3, pt.size = 0.1) * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10))
  SavePlot(od = singledir, filename = paste(sample.name, "filter.vln", sep = "_"), data = P.vlnplot)

  P1.scatter <- FeatureScatter(single.seurat, feature1 = "nUMI", feature2 = "percent.mt", pt.size = 0.5) + theme(legend.position = "none", axis.text = element_text(size = 10))
  P2.scatter <- FeatureScatter(single.seurat, feature1 = "nUMI", feature2 = "nGene", pt.size = 0.5) + theme(axis.text = element_text(size = 10))
  P.scatter <- P1.scatter + P2.scatter
  SavePlot(od = singledir, filename = paste(sample.name, "filter.scatter", sep = "_"), data = P.scatter)
  save(single.seurat, file = file.path(singledir, "single_qc.Rda"))
  saveRDS(single.seurat, file = file.path(singledir, "single_qc.Rds"))
  return(single.seurat)
}, mc.cores = length(object))
# })


names(singleList) <- sample
save(singleList,file = file.path(od, 'upload.Rdata'))

