library(getopt)
spec=matrix(c(
  'help','h',0,"logical",
  'matrix','m',1,"character",
  'species','s',1,'character',
  'od','o',1,'character'
),byrow=TRUE,ncol=4
)
opt=getopt(spec)
###help
print_usage<-function(spec=NULL){
  cat(getopt(spec,usage=TRUE));
  cat("Usage example: \n")
  cat("
		Description: Computational cell cycle, Draw related pictures
		Contact: liyu zheng <zhengly@biomarker>	
		Usage example:
			/share/nas2/genome/biosoft/R/3.6.1/bin/Rscript cell_cycle.R --matrix /path/single_seruat.Rds --species Mus_musculus --od outputs 
			/share/nas2/genome/biosoft/R/3.6.1/bin/Rscript cell_cycle.R -m /path/single_seruat.Rds -s Mus_musculus -o outputs
		
			
		Options:
		--help/-h	NULL		get this help
		--od/-o		character	output dir
    --matrix/-m  Single cell seurat data,*.Rds
            eg: single_seruat.Rds
    --species/-s  Corresponding species, human or mouse
		\n")
  
  q(status=1);
}
if(is.null(opt$m) && is.null(opt$s) && is.null(opt$o)) print_usage(spec)
if(opt$s!="Homo_sapiens" && opt$s!="Mus_musculus" && opt$s != "Rattus_norvegicus") print_usage(spec)
if(!exists(opt$o)) {dir.create(opt$o)}

library(Seurat)
library(tidyverse)
library(cowplot)
library(pheatmap)
seurat <- readRDS(opt$m)

###修改2021.12.13，路径问题，此处已经写死后期安装here包直接获取
file_path <- '/share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/cell_cycle/'
if(opt$s=="Homo_sapiens"){
  load(str_c(file_path,'Human_cell_cycle_gene.rda'))
  # load("/share/nas1/ranjr/pipeline/singleCell_trans/v1.0/cell_cycle/Human_cell_cycle_gene.rda")
}else{
  load(str_c(file_path,'Mouse_cell_cycle_gene.Rda'))
  # load("/share/nas1/ranjr/pipeline/singleCell_trans/v1.0/cell_cycle/Mouse_cell_cycle_gene.Rda")
}


counts <- seurat@assays$RNA@counts 
metadata <- seurat@meta.data
marrow <- CreateSeuratObject(counts = counts,meta.data = metadata)
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
marrow <- ScaleData(marrow, features = rownames(marrow))

marrow <- CellCycleScoring(marrow, s.features = s_genes, g2m.features = g2m_genes)
head(marrow[[]])
CC <- marrow@meta.data %>% rownames_to_column(var = "Cell") %>% select(Cell,sample,S.Score,G2M.Score,Phase)
write.table(CC,paste0(opt$o,"/cell_cycle.xls"),col.names = T,row.names = F,quote = F,sep = "\t")
sample_id <- unique(CC$sample)
CC_cluster <- marrow@meta.data %>% rownames_to_column(var = "Cell") %>% select(Cell,sample,Phase,seurat_clusters)

cluster_plot <- function(x){
  CC_cluster_test <-CC_cluster %>% filter(sample==x) %>% group_by(seurat_clusters,Phase) %>% summarise(con=n())
  p <- CC_cluster_test %>% ggplot(aes(x = seurat_clusters, y=con,fill=Phase)) + 
  geom_bar(stat = 'identity',position = position_dodge2(padding = 0,preserve = "single"))+
  theme_classic()+labs(x = '', y = 'number of cell', title = paste0(x,' Cell cycle information in each cluster'))
  scale_y_continuous(expand = c(0,0))
  ggsave(paste0(opt$o,"/",x,"_cell_cycle_cluster_barplot.pdf"),p,device ="pdf",units = "cm",
           dpi = 300,
           limitsize = TRUE)
  ggsave(paste0(opt$o,"/",x,"_cell_cycle_cluster_barplot.png"),p,device ="png",units = "cm",
           dpi = 300,
           limitsize = TRUE)
}
map(sample_id,cluster_plot)

# pie_plot <- function(x){
#     CC_A <- CC %>% filter(sample==x)
#     CC_A <- CC_A %>% group_by(sample,Phase) %>% summarise(con=n())
#     label_value <- paste('(', round(CC_A$con/sum(CC_A$con) * 100, 1), '%)', sep = '')
#     label <- paste(CC_A$Phase, label_value, sep = '')
#     p <-  CC_A %>% ggplot(mapping = aes(x = '', y = con, fill = Phase)) +
#       geom_bar(stat = 'identity', position = 'stack', width = 1)+coord_polar(theta = 'y')+
#       labs(x = '', y = '', title = paste0(x,' cell cycle relative abundance (%)')) +
#       theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank())+
#       scale_fill_discrete(labels = label)
#     ggsave(paste0(opt$o,"/",x,"cell_cycle_pie.pdf"),p,device ="pdf",units = "cm",
#            dpi = 300,
#            limitsize = TRUE)
#     ggsave(paste0(opt$o,"/",x,"cell_cycle_pie.png"),p,device ="png",units = "cm",
#            dpi = 300,
#            limitsize = TRUE)
# }
# map(sample_id,pie_plot)

# if(length(sample_id)==1){
#   label_value <- paste('(', round(CC_A$con/sum(CC_A$con) * 100, 1), '%)', sep = '')
#   label <- paste(CC_A$Phase, label_value, sep = '')
#   p <-  CC_A %>% ggplot(mapping = aes(x = '', y = con, fill = Phase)) + 
#     geom_bar(stat = 'identity', position = 'stack', width = 1)+coord_polar(theta = 'y')+
#     labs(x = '', y = '', title = 'Cell Cycle Relative abundance (%)') + 
#     theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank())+
#     scale_fill_discrete(labels = label)
#   ggsave(paste0(opt$o,"/cell_cycle_pie.pdf"),p,device ="pdf",units = "cm",
#          dpi = 300,
#          limitsize = TRUE)
#   ggsave(paste0(opt$o,"/cell_cycle_pie.png"),p,device ="png",units = "cm",
#          dpi = 300,
#          limitsize = TRUE)
# }else{
#   ####柱状堆砌�?
#   if(length(sample_id)==2){
#     p <- CC_A %>% ggplot(aes(x=sample,y=con,fill=Phase))+geom_bar(stat = "identity",position = "fill",width=0.6)+
#       xlab("")+ ylab("Cell Cycle Relative abundance (%)")+
#       theme_classic()+
#       coord_cartesian(ylim=c(0,1),xlim = c(1.2,(length(sample_id)-0.2)))+
#       scale_y_continuous(expand = c(0,0)) +
#       theme(aspect.ratio = 2/1)
#   }else{
#     p <- CC_A %>% ggplot(aes(x=sample,y=con,fill=Phase))+geom_bar(stat = "identity",position = "fill",width=0.6)+
#       xlab("")+ ylab("Cell Cycle Relative abundance (%)")+
#       theme_classic()+
#       coord_cartesian(ylim=c(0,1),xlim = c(1.2,(length(sample_id)-0.2)))+
#       scale_y_continuous(expand = c(0,0)) 
#   }
#   ggsave(paste0(opt$o,"/cell_cycle_barplot.pdf"),p,device ="pdf",units = "cm",
#          dpi = 300,
#          limitsize = TRUE)
#   ggsave(paste0(opt$o,"/cell_cycle_barplot.png"),p,device ="png",units = "cm",
#          dpi = 300,
#          limitsize = TRUE)
# }


###########pheatmap
##col_ann
col_ann <- CC %>% select(Cell,Phase)
col_ann <- col_ann[order(col_ann$Phase),] 
rownames(col_ann) <- col_ann$Cell
col_ann <- col_ann %>% select(-Cell)

###data
all_data <- data.frame(marrow@assays$RNA@scale.data,check.names = F)
s_p <- all_data %>% rownames_to_column(var = "gene") %>% 
  filter(gene %in% c(s_genes,g2m_genes)) %>% column_to_rownames(var = "gene")
##row_ann
row_ann <- data.frame(gene=rownames(s_p),
                      Class="")
row_ann <- row_ann %>% mutate(Class=case_when(
  gene %in% g2m_genes~"G2/M",
  gene %in% s_genes~"S"
)) 
row_ann <- row_ann[order(row_ann$Class),]
rownames(row_ann) <- row_ann$gene
row_ann <- row_ann %>% select(-gene)

s_p <- s_p[,rownames(col_ann)]
s_p <- s_p[rownames(row_ann),]

##设置颜色�?
#breaks
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

pheatmap(s_p,scale="row", color =c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),breaks=bk,
         show_rownames = F,show_colnames = F,
         cluster_cols = F,cluster_rows = F,
         annotation_col = col_ann,angle_col = "90",filename =paste0(opt$o,"/cell_cycle_heatmap.pdf"))
pheatmap(s_p,scale="row", color =c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),breaks=bk,
         show_rownames = F,show_colnames = F,
         cluster_cols = F,cluster_rows = F,
         annotation_col = col_ann,angle_col = "90",filename =paste0(opt$o,"/cell_cycle_heatmap.png"))
######pheatmap



