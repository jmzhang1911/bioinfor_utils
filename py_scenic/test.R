library(SCENIC)
library(tidyverse)
library(SCopeLoomR)
library(ComplexHeatmap)
library(RColorBrewer)
library(pheatmap)
library(circlize)
library(Cairo)

loom <- open_loom('./sample_SCENIC.loom')
glimpse(loom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC")
                               

n=t(scale(t( getAUC(regulonAUC[,] ))))
load('./meta_data.RData')
meta_data %>% rownames_to_column('barcodes') %>%
  filter(barcodes %in% colnames(n)) -> meta_df

getPalette_cell <- colorRampPalette(brewer.pal(8, "Set1"))
cell_ve <- unique(meta_df$cellType);ve_col <- getPalette_cell(length(cell_ve))
names(ve_col) <- cell_ve

top_anno_ref <- HeatmapAnnotation(
  celltype = meta_df$cellType, height = unit(0.1, 'npc'),
  col = list(celltype = ve_col)
)

ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), labels = month.name[1:10]))
col_fun = colorRamp2(c(-2, 0, 2), c("#2F4F4F", "white", "#FFA500"))



Heatmap(n, 
        row_names_gp = gpar(fontsize = 1),
        col = col_fun, 
        column_gap = unit(0, "mm"), 
        column_title_rot = 90,
        column_split  = as.factor(meta_df$cellType),
        column_title_gp  = gpar(fontsize = 8),
        show_column_names = F,
        show_row_dend = F,
        show_column_dend = F,
        top_annotation = top_anno_ref) ->  p


png(filename = 'tmp.png', width = 800, height = 800)
#draw(p, use_raster = TRUE, raster_quality = 30)
draw(p)
dev.off()
