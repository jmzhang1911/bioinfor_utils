library(Seurat)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

load('/Volumes/My_Passport/data_of_server/tmp_seob.RData')

tmp_seob@meta.data %>% rownames_to_column('barcodes') %>%
  filter(abs(S.Score)>0.15 | abs(G2M.Score) > 0.15) %>% pull(barcodes) -> prob
tmp_seob[, prob] -> seob

col_fun = colorRamp2(c(-2, 0, 2), c("#000080", "white", "#B22222"))
seob@assays$RNA@scale.data -> df
seob <- ScaleData(seob, assay = 'RNA')


Heatmap(df, 
        col = col_fun,
        show_row_names = F, 
        show_column_names = F,
        height = unit(15, 'cm'), 
        width = unit(15, 'cm'))


ggplot(tmp_seob@meta.data, aes(x = S.Score)) +
  geom_density()

ggplot(tmp_seob@meta.data, aes(x = G2M.Score)) +
  geom_density()
