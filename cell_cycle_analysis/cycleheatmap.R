#!/usr/bin/env Rscript
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
library(getopt)
spec=matrix(c(
  'sce'    , 'i', 1, "character", "single cell object",
  'od'     , 'o', 1, 'character', "output dir",
  'Species', 's', 1, 'character', 'species name',
  'help'   , 'h', 0, "logical"  , "print usage"
), byrow = TRUE, ncol = 5)

opt=getopt(spec)
###help
print_usage<-function(spec=NULL){
  cat(getopt(spec,usage=TRUE));
  cat("Usage example: \n")
  cat("
                Description: Computational cell cycle, Draw related pictures
                Contact: liyu zheng <zhengly@biomarker>

                Options:
                --help/-h       NULL            get this help
                --od/-o         character       output dir
    --sce/-i  Single cell seurat data,*.Rds
            eg: single_seruat.Rds
    --Species/-s  Corresponding species, human or mouse
                \n")

  q(status=1);
}

if(is.null(opt$sce)) {print_usage(spec)}
if(is.null(opt$od)) {opt$od <- "./"}
if(!exists(opt$od)) {dir.create(opt$od)}

library(Seurat)
library(tidyverse)
library(cowplot)
library(pheatmap)
seurat.obj <- readRDS(opt$sce)

if(all(c("S.Score", "G2M.Score", "Phase") %in% (seurat.obj@meta.data %>% colnames))){


prefix <- unlist(strsplit(basename(opt$sce), split = "[.]"))[1]
od <- opt$od
my.Species <- opt$Species
if(my.Species == "human" || my.Species == "Homo_sapiens" || my.Species == "Human"){
  load(file.path(script.basename, "Human_cell_cycle_gene.Rda"))
}else if(my.Species == "mouse" || my.Species == "Mus_musculus" || my.Species == "Mouse"){
  load(file.path(script.basename, "Mouse_cell_cycle_gene.Rda"))
}else{
  q()
}
 
p <- DimPlot(seurat.obj, reduction = "umap", group.by = "Phase")
P.umap <- UMAPPlot(seurat.obj, group.by = "Phase", label = FALSE) + labs(title = "")
P.tsne <- TSNEPlot(seurat.obj, group.by = "Phase", label = FALSE) + labs(title = "")
ggsave(file = file.path(od, paste(prefix, "cellCycle_umap.png", sep = ".")), P.umap, width = 6, height = 4, scale = 1.3)
ggsave(file = file.path(od, paste(prefix, "cellCycle_umap.pdf", sep = ".")), P.umap, width = 6, height = 4, scale = 1.3)
ggsave(file = file.path(od, paste(prefix, "cellCycle_tsne.png", sep = ".")), P.tsne, width = 6, height = 4, scale = 1.3)
ggsave(file = file.path(od, paste(prefix, "cellCycle_tsne.pdf", sep = ".")), P.tsne, width = 6, height = 4, scale = 1.3)
 
CC <- seurat.obj@meta.data %>% rownames_to_column(var = "Cell") %>% select(Cell, sample, S.Score, G2M.Score, Phase)
write.table(CC, file.path(od, paste(prefix, "cell_cycle.xls", sep = ".")), col.names = T, row.names = F, quote = F, sep = "\t")
sample_id <- unique(CC$sample)
CC_cluster <- seurat.obj@meta.data %>% rownames_to_column(var = "Cell") %>% select(Cell, sample, Phase, seurat_clusters)

###########pheatmap
##col_ann
col_ann <- CC %>% select(Cell, Phase)
col_ann <- col_ann[order(col_ann$Phase),] 
rownames(col_ann) <- col_ann$Cell
col_ann <- col_ann %>% select(-Cell)

###data
all_data <- data.frame(seurat.obj@assays$RNA@data, check.names = F)
s_p <- all_data %>% rownames_to_column(var = "gene") %>% filter(gene %in% c(s_genes, g2m_genes)) %>% column_to_rownames(var = "gene")
##row_ann
row_ann <- data.frame(gene=rownames(s_p), Class = "")
row_ann <- row_ann %>% mutate(Class = case_when(gene %in% g2m_genes~"G2/M", gene %in% s_genes~"S")) 
row_ann <- row_ann[order(row_ann$Class),]
rownames(row_ann) <- row_ann$gene
row_ann <- row_ann %>% select(-gene)

s_p <- s_p[,rownames(col_ann)]
s_p <- s_p[rownames(row_ann),]

##设置颜色条
#breaks
bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))

pheatmap(s_p, scale = "row", color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2, 2, 1), breaks = bk,
         show_rownames = F, show_colnames = F,
         cluster_cols = F, cluster_rows = F,
         annotation_col = col_ann,angle_col = "90",filename = file.path(od, paste(prefix, "cell_cycle_heatmap.pdf", sep = ".")))
pheatmap(s_p, scale = "row", color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2, 2, 1), breaks = bk,
         show_rownames = F, show_colnames = F,
         cluster_cols = F, cluster_rows = F,
         annotation_col = col_ann, angle_col = "90", filename = file.path(od, paste(prefix, "cell_cycle_heatmap.png", sep = ".")))


}

Sys.sleep(150)
