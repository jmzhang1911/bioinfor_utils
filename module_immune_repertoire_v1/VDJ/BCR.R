#!/usr/bin/Rscript
# parameter
suppressMessages(library(getopt))
spec <- matrix(c(
  'input'          , 'i', 2, 'character', 'input dir',
  'od'             , 'o', 2, 'character', 'output dir',
  'type'           , 't', 2, 'character', 'chain of TCR',
  'species'        , 's', 2, 'character', 'species',
  'help'           , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:BCR.R
Version: Version v1.0
Contact: zhengly <zhengly@biomarker.com.cn> 
Program Date:   2020.12.25 
Description: this program is used to BCR ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$input) | is.null(opt$od) | is.null(opt$type)) {PrintUsage()}
if(!dir.exists(opt$od)) {dir.create(opt$od)}
od=opt$od
input=opt$input
type=opt$type
species=opt$species

##load packages
#.libPaths("/share/nas1/zhengly/softwares/R/")
.libPaths("/share/nas2/genome/biosoft/R/4.0.3/lib64/R/library.bak_20210428/")
suppressMessages(library(tidyverse))
suppressMessages(library(immunarch))
suppressMessages(library(ggsci))
suppressMessages(library(venn))
suppressMessages(library(patchwork))
suppressMessages(library(circlize))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggforce))

time1 <- proc.time()
dir.create(paste0(od,"/in"))
files <- list.files(input)
for (i in files) {
  data <- read_csv(paste0(input,"/",i))
  TR <- data %>% filter(chain==type)
  if(nrow(TR)>0){
  write.table(TR,paste0(od,"/in/",i),col.names = T,row.names = F,quote = F,sep = ",")
  }
}

##load data
immdata <- repLoad(paste0(od,"/in"),.mode = "single")

##function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}

##CDR3 lens
len_nt <- repExplore(immdata$data,.col = "aa",.method = "len")

len_nt_p <- ggplot(len_nt) + geom_col(aes(x = Length, y = Count,fill=Sample),
  position = position_dodge2(padding = 0, preserve = "single"))+
  theme_bw()+theme(panel.grid =element_blank())+scale_fill_npg()+
  labs(x="CDR3(aa) length",y="Clonotypes",fill = "Sample",title = paste0("Distribution of ",type," CDR3 polypeptide lengths"))+
  theme(plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=seq(0, max(len_nt$Length), 2))
SavePlot(od,paste0(type,"_CDR3_lens"),len_nt_p)
print("CDR3 lens done..")

##relative abundance
clon_rel <- repClonality(immdata$data, "homeo")
clon_rel_p <- clon_rel %>% vis()+scale_fill_npg()
SavePlot(od,paste0(type,"_Clonotype_Relative_abundance"),clon_rel_p)
print("relative abundance done..")

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
write.table(all_e,paste0(od,"/",type,"_diversity_measure.xls"),col.names = T,row.names = F,quote = F,sep="\t")
print("diversity measure done..")

##trackClonotypes
if(n_distinct(immdata$meta$Sample)>1){
track <- trackClonotypes(immdata$data, list(1, 5), .col = "aa") %>% vis()
SavePlot(od,paste0(type,"_trackClonotypes"),track)
print("trackClonotypes done..")
}

##V/J Gene
vizGene <- function(immdata=immdata,Gene="ighV"){
  if(str_sub(Gene,4,4)=="V"){
    gene <- geneUsage(immdata[["data"]], "hs.ighv")
  }else if(str_sub(Gene,4,4)=="J"){
    gene <- geneUsage(immdata[["data"]], "hs.ighj")
  }
  trav <- gene[grepl(x = gene$Names,pattern = Gene),]
  trav[is.na(trav)] <- 0
  if(n_distinct(immdata$meta$Sample)>1){
    trav <- trav %>% group_by(Names) %>% summarise(across(immdata$meta$Sample,sum))
    for(i in immdata$meta$Sample){
      trav[i] <- trav[i]/sum(trav[i])
    }
    trav <- trav %>% pivot_longer(immdata$meta$Sample,names_to = "Sample",values_to = "Count")
  }else{
    trav <- trav %>% group_by(Names) %>% summarise(Clones=sum(Clones))
    trav$Count <- trav$Clones/sum(trav$Clones)
    trav <- mutate(trav,Sample=immdata$meta$Sample)
  }
  p <- trav %>% ggplot() + geom_col(aes(x=Names,y=Count,fill=Sample),
                                    position =position_dodge2(padding = 0, preserve = "single"))+
    theme_bw()+theme(panel.grid =element_blank())+scale_fill_npg()+
    labs(x=paste(str_sub(Gene,1,3),str_sub(Gene,4,4),"Gene",collapse = ""),y="Percentage(%)",fill = "Sample",title = "")+
    theme(axis.text.x=element_text(angle=90, hjust=1))
  return(p)
}

if(type=="IGH"){
  TRAV=vizGene(immdata,"IGHV")
  SavePlot(od,paste0(type,"_IGHV_gene"),TRAV,width = 12,height = 8)
  TRAJ=vizGene(immdata,"IGHJ")
  SavePlot(od,paste0(type,"_IGHJ_gene"),TRAJ,width = 12,height = 8)
}else if(type=="IGK"){
  TRBV=vizGene(immdata,"IGKV")
  SavePlot(od,paste0(type,"_IGKV_gene"),TRBV,width = 12,height = 8)
  TRBJ=vizGene(immdata,"IGKJ")
  SavePlot(od,paste0(type,"_IGKJ_gene"),TRBJ,width = 12,height = 8)
}else if(type=="IGL"){
  TRBV=vizGene(immdata,"IGLV")
  SavePlot(od,paste0(type,"_IGLV_gene"),TRBV,width = 12,height = 8)
  TRBJ=vizGene(immdata,"IGLJ")
  SavePlot(od,paste0(type,"_IGLJ_gene"),TRBJ,width = 12,height = 8)
}
print("V/J Gene done..")
##V-J
for (i in 1:length(immdata$meta$Sample)) {
  V_J <- immdata$data[[i]] %>%
    select(V.name,J.name,Clones)%>% 
    group_by(V.name,J.name) %>% summarise(count=sum(Clones))
  V_J_p  <-V_J %>% ggplot(aes(x=V.name,y=J.name))+
    geom_tile(alpha=0.5,aes(fill=count))+
    theme(axis.text.x=element_text(angle=90,vjust=0.5))+    
    labs(x="V Gene",y="J Gene",title=" ",fill="Clonotypes")+   
    theme(panel.grid = element_blank())+            
    scale_fill_gradient2(midpoint =mean(c(min(V_J$count),max(V_J$count))),
                         low = "#6181BD4E",mid = "#64A10E4E",high = "#FDA1004E")
  SavePlot(od,paste0(type,"_V_J_pair_",immdata$meta$Sample[i]),V_J_p,height = 8,width = 12)
}
print("V/J pair done..")
##crioc
circos_plot <- function(mat,od=od,name="circos_plot.png"){
    if (grepl("\\.pdf$",name)){
       pdf(paste0(od,"/",name))
    } else if (grepl("\\.png$",name)) {
       png(paste0(od,"/",name), 
           width     = 3.25,
           height    = 3.25,
           units     = "in",
           res       = 1200,
           pointsize = 4)
    } 
    mat <- as.matrix(mat)
    col_sum = apply(mat, 2, sum)
    row_sum = apply(mat, 1, sum)
    mat <- mat[order(row_sum), order(col_sum)]
    rn <- rownames(mat)
    cn <- colnames(mat)
    maxrn <- max(nchar(rn))
    maxcn <- max(nchar(cn))
    for(i in seq_len(length(rn))) {
      rn[i] <- paste(rn[i], paste(rep(" ", maxrn - nchar(rn[i])), collapse = ''))
    }
    for(i in seq_len(length(cn))) {
      cn[i] <- paste(cn[i], paste(rep(" ", maxcn - nchar(cn[i])), collapse = ''))
    }
    rownames(mat) <- rn
    colnames(mat) <- cn
    circos.par(gap.degree = c(rep(1, nrow(mat)-1), 10, rep(1, ncol(mat)-1), 15), start.degree = 5)
    rcols <- rep(brewer.pal(12, "Paired"), nrow(mat)/12 + 1)[1:nrow(mat)]
    ccols <- rep(brewer.pal(12, "Paired"), ncol(mat)/12 + 1)[1:ncol(mat)]
    names(rcols) <- sort(rownames(mat))
    names(ccols) <- sort(colnames(mat))
    chordDiagram(mat, annotationTrack = "grid",
             grid.col = c(rcols, ccols),
             preAllocateTracks = list(track.height = 0.2), transparency = 0.5)
    circos.trackPlotRegion(track.index = 1, bg.border = NA,
                       panel.fun = function(x, y) {
                         sector.name = get.cell.meta.data("sector.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                       }
    )
    circos.clear()
    dev.off()
}
for(i in 1:length(immdata$meta$Sample)){
  test <- immdata$data[[i]] %>% select(V.name,J.name,Clones)
  test <- test %>% group_by(V.name,J.name) %>% summarise(Clones=sum(Clones))
  test_ <- test %>% pivot_wider(names_from = V.name,values_from=Clones)
  test_ <- test_ %>% column_to_rownames(var = "J.name")
  test_[is.na(test_)] <- 0
  circos_plot(test_,od,paste0(type,"_Circos_",immdata$meta$Sample[i],".pdf"))
  circos_plot(test_,od,paste0(type,"_Circos_",immdata$meta$Sample[i],".png"))
}
print("Circos done..")
##V-D-J
for(i in 1:length(immdata$meta$Sample)){
    test <- immdata$data[[i]] %>%
      select(V.name,D.name,J.name,Clones) %>% 
      filter(!D.name %in% c("NA","none","None"))
    test <-na.omit(test)
    if(nrow(test)>5){
        test <- test %>% group_by(V.name,D.name,J.name) %>% summarise(count=sum(Clones))
        p <- test%>%
            gather_set_data(c(2,1,3)) %>%
            ggplot(aes(x, id = id, split = y, value = 1))  +
            geom_parallel_sets(aes(fill = J.name), show.legend = FALSE, alpha = 0.3) +
            geom_parallel_sets_axes(axis.width = 0.1, color = "lightgrey", fill = "white") +
            geom_parallel_sets_labels(angle = 0) +
            theme_no_axes()
            SavePlot(od,paste0(type,"_V-D-J_",immdata$meta$Sample[i]),p,width = 12,height = 8)
    }
}
print("V-D-J done")
###V/J in CDR3
for(i in 1:length(immdata$meta$Sample)){
  spect_v <- spectratype(immdata$data[[i]], .quant = "count", .col = "aa+v") 
  if(nrow(spect_v)>15 && sum(!duplicated(spect_v$Gene))>=12){
    spect_v <- spect_v %>% vis()
    SavePlot(od,paste0(type,"_CDR3_V_",immdata$meta$Sample[i]),spect_v,width = 12,height = 8)
  }
}
print("V/J CDR3 done..")

###motif
for(i in 1:length(immdata$meta$Sample)){
  les <- min(str_length(immdata$data[[i]]$CDR3.aa))
  kmers <- getKmers(immdata$data[[i]], les,.coding = T,.col="aa")
  p <- kmer_profile(kmers, "self") %>% vis(.plot = "seq")
  SavePlot(od,paste0(type,"_motif_",immdata$meta$Sample[i]),p)
}
print("motif done..")



time2 <- proc.time()
run.time <- time2 - time1
print(paste0('Elapsed time: ',run.time[3][[1]],"s"))

















