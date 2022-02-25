#!/usr/bin/Rscript
# parameter
suppressMessages(library(getopt))
spec <- matrix(c(
  'input'          , 'i', 2, 'character', 'input dir',
  'od'             , 'o', 2, 'character', 'output dir',
  'rds'            , 'r', 2, 'character', 'seurat.rda',
  'type'           , 't', 2, 'character', 'B or T',
  'help'           , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:scRepertoire.R
Version: Version v1.0
Contact: zhengly <zhengly@biomarker.com.cn> 
Program Date:   2021.03.23 
Description: use scRepertoire analysis sc+vdj ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$input) | is.null(opt$od) | is.null(opt$type) |is.null(opt$rds)) {PrintUsage()}
if(!dir.exists(opt$od)) {dir.create(opt$od)}
od=opt$od
path=opt$input
type=opt$type

time1 <- proc.time()
.libPaths("/share/nas2/genome/biosoft/R/4.0.3/lib64/R/library.bak_20210428/")
suppressMessages(library(scRepertoire))
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}
###数据读取
file_name <- list.files(path)
contig_lists <-c()
for(i in seq_along(file_name)){
  contig_lists[[i]] <- read_csv(paste0(path,"/",file_name[i]))
}
###多样本的结果合并 这边分TCR和BCR
file_name=sub(file_name,pattern = ".csv",replacement = "")
sample_id=t(data.frame(str_split(file_name,pattern = "-")))
if(type=="B"){
  combined <- combineBCR(contig_lists, 
                         samples = unname(sample_id[,1]), 
                         ID = unname(sample_id[,2]))
  
}
if(type=="T"){
  combined <- combineTCR(contig_lists,
                         samples = unname(sample_id[,1]),
                         ID = unname(sample_id[,2]), 
                         cells ="T-AB")
}

###读取单细胞的seurat分析结果  这边多样本使用多样本的结果去合并
seurat <- get(load(opt$rds))
cluster_u <- UMAPPlot(seurat)
cluster_t <- TSNEPlot(seurat)
###barcode信息处理
for (i in seq_along(combined)) {
  combined[[i]] <- stripBarcode(combined[[i]] , column = 1, connector = "_", num_connects = 3)
}
if(length(combined)>1){
  sample_name <- seurat@meta.data %>% rownames_to_column(var = "Cell") %>% select(Cell,sample)
  sample_name <- sample_name %>% rowwise() %>% mutate(id=str_split(Cell,"_")[[1]][2]) %>% select(sample,id) %>% unique()
  sample_name$sample <- sub(sample_name$sample,pattern = "-sc",replacement = paste0("_",unique(unname(sample_id[,2]))))
  for (i in sample_name$sample) {
    id=which(sample_name$sample==i)
    combined[[i]]$barcode <-paste0(combined[[i]]$barcode,"_",sample_name$id[id]) 
  }
}
###VDJ与seurat数据进行合并
seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene+nt", groupBy = "sample",
                            cloneTypes=c(None = 0,Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
#cloneType_p <- DimPlot(seurat, group.by = "cloneType")
#SavePlot(od=od,filename = "seurat_cloneType",data = cloneType_p)
cloneType_ump <- UMAPPlot(seurat, group.by = "cloneType") + labs(color = "Clonotype group",title = "")
cloneType_ump_a <- cluster_u + cloneType_ump
SavePlot(od=od,filename = paste0(type,"CR_cloneType_ump"),data = cloneType_ump_a, width = 12)
cloneType_tsne <- TSNEPlot(seurat, group.by = "cloneType") + labs(color = "Clonotype group",title = "")
cloneType_tsne_a <- cluster_t + cloneType_tsne
SavePlot(od=od,filename = paste0(type,"CR_cloneType_tsne"),data = cloneType_tsne_a, width = 12)

cluster_p<- occupiedscRepertoire(seurat, x.axis = "cluster")+labs(fill = "Clonotype group")+ylab("Cells")
SavePlot(od=od,filename = paste0(type,"CR_cluster_cloneType"),data = cluster_p)

combined2 <- expression2List(seurat, group = "cluster")
#cluster_Diversity<- clonalDiversity(combined2, cloneCall = "aa")
#SavePlot(od=od,filename = "cluster_Diversity",data = cluster_Diversity)
#clonalDiversity(combined2,exportTable = T)

#cluster_clonal <- clonalHomeostasis(combined2, cloneCall = "aa")
#SavePlot(od=od,filename = "cluster_clonal",data = cluster_clonal)

cluster_Overlap <- clonalOverlap(combined2, cloneCall="aa", method="overlap")
SavePlot(od=od,filename =paste0(type,"CR_cluster_Overlap"),data = cluster_Overlap)

a <- seurat@meta.data %>% group_by(sample) %>% arrange(desc(Frequency)) %>% summarise(CTaa=unique(CTaa)) %>% slice(1:5)
a <- a  %>% na_if("NA_NA") %>% na.omit()
b <- seurat@meta.data
b$clones <- ""
for (i in seq_len(nrow(a))) {
  b <- b %>% mutate(clones=case_when(
    sample==a$sample[i] & CTaa==a$CTaa[i] ~ paste0(a$sample[i],":",a$CTaa[i]),
    TRUE~clones
  ))
}
b$clones <- na_if(b$clones,"")
seurat@meta.data <- b
clones_ump <- UMAPPlot(seurat, group.by = "clones") + labs(color = "The largest number of clones",title = "")
#SavePlot(od=od,filename = paste0(type,"CR_Clones_ump"),data = clones_ump, width = 12)
clones_tsne <- TSNEPlot(seurat, group.by = "clones") + labs(color = "The largest number of clones",title = "")
#SavePlot(od=od,filename = paste0(type,"CR_Clones_tsne"),data = clones_tsne, width = 12)

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('Elapsed time: ',run.time[3][[1]],"s"))











