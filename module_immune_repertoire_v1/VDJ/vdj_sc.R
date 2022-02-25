#!/usr/bin/Rscript
# parameter
suppressMessages(library(getopt))
spec <- matrix(c(
  'Tpath'          , 't', 2, 'character', 'T input dir',
  'od'             , 'o', 2, 'character', 'output dir',
  'Bpath'        , 'b', 2, 'character', 'B input dir',
  'help'           , 'h', 0, 'logical'  , 'print usage',
  'rda'           , 'r', 2, 'character'  , 'single cell data'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:vdj_sc.R
Version: Version v1.0
Contact: zhengly <zhengly@biomarker.com.cn> 
Program Date:   2021.01.10 
Description: this program is used to vdj+sc ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help)) {PrintUsage()}
if((is.null(opt$Tpath) & is.null(opt$Bpath)) | is.null(opt$od) | is.null(opt$rda)) {PrintUsage()}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

time1 <- proc.time()
##load packages
suppressMessages(library(Seurat,lib.loc ="/share/nas1/zhengly/softwares/R/" ))
suppressMessages(library(tidyverse))
##function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}
od <- opt$od

load(opt$rda)

add_vdj <- function(data=data,type=type,input=input){
  tcr <- read.csv(input)
  tcr <- tcr%>% group_by(barcode) %>% 
    summarise(chain=paste(sort(unique(chain)),collapse = ";"),
              v_gene=paste(sort(unique(v_gene)),collapse = ";"),
              d_gene=paste(sort(unique(d_gene)),collapse = ";"),
              j_gene=paste(sort(unique(j_gene)),collapse = ";"),
              clonotype_id=paste(sort(unique(raw_clonotype_id)),collapse = ";")) %>% 
    mutate(chain=case_when(
      chain %in% "Mulit"~"Multi",
      TRUE~chain
    ))
  clono_file <- gsub(input,pattern = "filtered_contig_annotations.csv",replacement = "clonotypes.csv")
  clono <- read.csv(clono_file)
  tcr <- merge(tcr, clono,by = "clonotype_id")
  #Reorder so barcodes are first column and set them as rownames.
  rownames(tcr) <- tcr[,2]
  tcr[,2] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=data, metadata=tcr)
  return(clono_seurat)
}
if(!is.null(opt$Tpath)){
  single.seurat <- add_vdj(data = single.seurat,type = "t",input = opt$Tpath)
}
if(!is.null(opt$Bpath)){
  single.seurat <- add_vdj(data = single.seurat,type = "b",input = opt$Bpath)
}
#table(!is.na(single.seurat$t_clonotype_id),!is.na(single.seurat$b_clonotype_id))
if(!is.null(opt$Tpath)&!is.null(opt$Bpath)){
  single.seurat<- subset(single.seurat, cells = colnames(single.seurat)[!(!is.na(single.seurat$t_clonotype_id) & !is.na(single.seurat$b_clonotype_id))])
}
#table(!is.na(single.seurat$t_clonotype_id),!is.na(single.seurat$b_clonotype_id))

p1 <- DimPlot(single.seurat,reduction = "tsne", label = TRUE)
p4 <- DimPlot(single.seurat,reduction = "umap", label = TRUE)
if(!is.null(opt$Tpath)){
  p2 <-DimPlot(single.seurat,reduction = "tsne",group.by = "t_chain")
  p12<-p1+p2
  SavePlot(od,"t_chain_tsne",p12)
  p5 <-DimPlot(single.seurat,reduction = "umap",group.by = "t_chain")
  p45<-p4+p5
  SavePlot(od,"t_chain_umap",p45) 
}

if(!is.null(opt$Bpath)){
  p3 <-DimPlot(single.seurat,reduction = "tsne",group.by = "b_chain")
  p13<-p1+p3
  SavePlot(od,"b_chain_tsne",p13)
  p6 <-DimPlot(single.seurat,reduction = "umap",group.by = "b_chain")
  p46<-p4+p6
  SavePlot(od,"b_chain_umap",p46)
}
save(single.seurat, file = file.path(od, "vdj_sc.Rda"))
saveRDS(single.seurat, file = file.path(od, "vdj_sc.Rds"))

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('Elapsed time: ',run.time[3][[1]],"s"))

