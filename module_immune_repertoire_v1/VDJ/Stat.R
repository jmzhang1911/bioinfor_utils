#!/usr/bin/Rscript
# parameter
suppressMessages(library(getopt))
spec <- matrix(c(
  'input'          , 'i', 2, 'character', 'input dir',
  'od'             , 'o', 2, 'character', 'output dir',
  'type'           , 't', 2, 'character', 'type',
  'help'           , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:Stat.R
Version: Version v1.0
Contact: zhengly <zhengly@biomarker.com.cn> 
Program Date:   2021.01.08 
Description: this program is used to statistical information;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$input) | is.null(opt$od)) {PrintUsage()}
if(!dir.exists(opt$od)) {dir.create(opt$od)}
od=opt$od
input=opt$input
type=opt$type
.libPaths("/share/nas1/zhengly/softwares/R/")
suppressMessages(library(tidyverse))
suppressMessages(library(immunarch))
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(type,filename, "png", sep = ".")
  file.pdf <- paste(type,filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}

immdata <- repLoad(input,.mode = "single")
cdr3 <- c()
for (i in 1:length(immdata$meta$Sample)) {
  aa <- length(unique(immdata$data[[i]]$CDR3.aa))
  v <- length(unique(immdata$data[[i]]$V.name))
  j <- length(unique(immdata$data[[i]]$J.name))
  vj <- length(unique(paste0(immdata$data[[i]]$V.name,immdata$data[[i]]$J.name)))
  mid <- immdata$data[[i]] %>% filter(!is.na(D.name)) %>% filter(!D.name %in% c("NA","none","None"))
  d <- length(unique(mid$D.name))
  vdj <-  length(unique(paste0(mid$V.name,mid$J.name)))
  cdr3 <- rbind(cdr3,c(immdata$meta$Source[i],immdata$meta$Chain[i],aa,v,d,j,vj,vdj))
}
colnames(cdr3) <- c("Sample","Chain","CDR3","V_Gene_variety","D_Gene_variety","J_Gene_variety","V-J_pair","V-D-J_pair")
write.table(cdr3,paste0(od,"/",type,".All_Samples_stat.xls"),col.names = T,row.names = F,quote = F,sep="\t")

rm(immdata)

immdata <- repLoad(input,.mode = "paired")
##Sample correlation
if(n_distinct(immdata$meta$Sample)>2){
imm_over_p <- repOverlap(immdata$data,.method = "overlap",.col = "aa", .verbose = F) %>% vis(.text.size=2)+
              labs(x="",y="")+theme(plot.title = element_text(hjust = 0.5))
#imm_over_p$layers[[2]]$aes_params$size <- 3
SavePlot(od,"All_Sample_overlap",imm_over_p)
print("Sample correlation done..")
}


##diversity measure
chao1_ <- repDiversity(immdata$data, "chao1",.col = "aa") %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% 
  select(Sample,Estimator) %>% 
  rename(Chao1=Estimator)

simp_ <- repDiversity(immdata$data, "gini.simp",.col = "aa") %>% 
  rename(Gini_Simpson=Value)
  
d50_ <- repDiversity(immdata$data, "d50",.col = "aa") %>% 
  as.data.frame() %>% 
  rownames_to_column(var="Sample") %>% 
  select(Sample,Clones) %>% 
  rename(D50=Clones) 

gini_ <- repDiversity(immdata$data, "gini",.col = "aa")%>% 
  as.data.frame() %>% 
  rename(Gini=V1) %>% 
  rownames_to_column(var="Sample")

all_e <- inner_join(chao1_,gini_,"Sample") %>% 
  inner_join(.,simp_,"Sample") %>% 
  inner_join(.,d50_,"Sample")
write.table(all_e,paste0(od,"/",type,".All_Sample_diversity_measure.xls"),col.names = T,row.names = F,quote = F,sep="\t")  
print("diversity measure done..")

###rarefaction
step <- if(min(unlist(lapply(immdata$data, function(x) sum(x$Clones))))<50) 1 else  NA
rare <- repDiversity(immdata$data, .method = "raref",.col = "aa", .step=step) %>% vis()
SavePlot(od,"All_Sample_rarefaction",rare)
print("rarefaction done..")

##Rank-abundance
rank_p <- repExplore(immdata$data, .method = "count") %>% vis()
SavePlot(od,"All_Sample_Rank-abundance",rank_p)
print("Rank-abundance done..")

#V_cor;cluster
vGene <- function(immdata =immdata,Gene="TRAV",od=od){
  imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T) 
  imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F) %>%
    vis(., .leg.title = "Cor", .text.size = 2.5)+
    labs(x="",y="")
  SavePlot(od,"All_Sample_VGene_cor",imm_gu_cor)
  imm_gu_pca <- geneUsageAnalysis(imm_gu, "pca",.verbose = F) %>% vis()
  SavePlot(od,"All_Sample_VGene_pca",imm_gu_pca)
  imm_gu_clust <- geneUsageAnalysis(imm_gu, "js+hclust",.verbose = F) %>% vis() 
  SavePlot(od,"All_Sample_VGene_clust",imm_gu_clust)
}


if(n_distinct(immdata$meta$Sample)>2){
  vGene(immdata,"TRAV",od)
}










