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
Program Date:   2021.01.22 
Description: this program is used to vdj+sc ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help)) {PrintUsage()}
if((is.null(opt$Tpath) & is.null(opt$Bpath)) | is.null(opt$od) | is.null(opt$rda)) {PrintUsage()}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

time1 <- proc.time()
##load packages
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
##function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}
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

od <- opt$od
single.seurat <- readRDS(opt$rda)
get_names <- function(x){
  a=basename(strsplit(x,"/outs/filtered_contig_annotations.csv")[[1]][1])
  a=unlist(strsplit(a,"-"))[1]
  return(a)
}
if(!is.null(opt$Tpath)){
  Tpath <- unlist(strsplit(opt$Tpath,","))
  samples <- unlist(lapply(Tpath,get_names))
}
if(!is.null(opt$Bpath)){
  Bpath <- unlist(strsplit(opt$Bpath,","))
  samples <- unlist(lapply(Bpath,get_names))
}

for(i in samples){
  samples.seurat <- subset(single.seurat,subset = orig.ident == paste0(i,"-sc"))
  rownam <- rownames(samples.seurat@meta.data)
  rownames(samples.seurat@meta.data) <- unlist(map(rownames(samples.seurat@meta.data),function(x){
    a=unlist(strsplit(x,"_"))[1]
    return(a)
  }))
  if(!is.null(opt$Tpath)){
    t_path <- Tpath[grep(Tpath,pattern =paste0(i,"-t/outs/filtered_contig_annotations.csv"))]
    samples.seurat <- add_vdj(data = samples.seurat,type = "t",input = t_path)
  }
  if(!is.null(opt$Bpath)){
    b_path <- Bpath[grep(Bpath,pattern =paste0(i,"-b/outs/filtered_contig_annotations.csv"))]
    samples.seurat <- add_vdj(data = samples.seurat,type = "b",input = b_path)
  }
  if(!is.null(opt$Tpath)&!is.null(opt$Bpath)){
    samples.seurat<- subset(samples.seurat, cells = colnames(samples.seurat)[!(!is.na(samples.seurat$t_clonotype_id) & !is.na(samples.seurat$b_clonotype_id))])
  }
  rownames(samples.seurat@meta.data) <- rownam
  p1 <- DimPlot(samples.seurat,reduction = "tsne", label = TRUE)
  p4 <- DimPlot(samples.seurat,reduction = "umap", label = TRUE)
  if(!is.null(opt$Tpath)){
    p2 <-DimPlot(samples.seurat,reduction = "tsne",group.by = "t_chain")
    p12<-p1+p2
    SavePlot(od,paste0(i,"-t_chain_tsne"),p12)
    p5 <-DimPlot(samples.seurat,reduction = "umap",group.by = "t_chain")
    p45<-p4+p5
    SavePlot(od,paste0(i,"-t_chain_umap"),p45) 
  }
  
  if(!is.null(opt$Bpath)){
    p3 <-DimPlot(samples.seurat,reduction = "tsne",group.by = "b_chain")
    p13<-p1+p3
    SavePlot(od,paste0(i,"-b_chain_tsne"),p13)
    p6 <-DimPlot(samples.seurat,reduction = "umap",group.by = "b_chain")
    p46<-p4+p6
    SavePlot(od,paste0(i,"-b_chain_umap"),p46)
  }
  save(samples.seurat, file = file.path(od, paste0(i,"-sc_vdj_sc.Rda")))
  saveRDS(samples.seurat, file = file.path(od, paste0(i,"-sc_vdj_sc.Rds")))
}
time2 <- proc.time()
run.time <- time2 - time1
print(paste0('Elapsed time: ',run.time[3][[1]],"s"))





























