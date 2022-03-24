library(SCENIC)
library(SCopeLoomR)
library(ComplexHeatmap)
library(pheatmap)
library(Cairo)

loom <- open_loom('./sample_SCENIC.loom')
glimpse(loom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC")
                               

n=t(scale(t( getAUC(regulonAUC[,] ))))

Heatmap(n, show_column_names = F)
