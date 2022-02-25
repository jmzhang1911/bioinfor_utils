MyBasicAnalysisSingleCell <- function(SeuratObj_, res_=0.8){
    SeuratObj_ <- NormalizeData(SeuratObj_, normalization.method = "LogNormalize")
    SeuratObj_ <- FindVariableFeatures(SeuratObj_, selection.method = "vst", nfeatures = 2000)
    SeuratObj_ <- ScaleData(SeuratObj_)
    
    SeuratObj_ <- RunPCA(SeuratObj_, features = VariableFeatures(SeuratObj_))
    SeuratObj_ <- RunUMAP(SeuratObj_, dims = 1:30)
    SeuratObj_ <- RunTSNE(SeuratObj_, dims = 1:30)
    
    SeuratObj_ <- FindNeighbors(SeuratObj_, dims = 1:30)
    SeuratObj_ <- FindClusters(SeuratObj_, resolution = res_)

    return(SeuratObj_)
}

MyRunDoubletFinderGetPK <- function(sweep.res.list_,outputdir_='./',title='GetPK'){
    sweep.stats <- summarizeSweep(sweep.res.list_, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    bcmvn %>% slice_max(BCmetric, n = 1) %>% pull(pK) -> x_intercept
    pk_good <- as.numeric(as.character(x_intercept))
    bcmvn %>%  
    ggplot(aes(x = pK, y = BCmetric)) +
        geom_point() +
        geom_line(aes(group = 1)) + 
        geom_vline(xintercept = x_intercept, color = 'red', linetype = 'dashed') +
        theme_classic() +
        theme(text = element_text(size = 18, face = 'bold'),
              axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) -> p
    MyPlotsSave(p, filename = title,output = outputdir_, width = 8, height = 6)
    
    return(pk_good)
}


MyRunDoubletFinder <- function(seuratObj_,cluster_colname_='seurat_clusters',pk_good,output_='./',title_='run_DoubletFinder_'){
    if(cluster_colname_ != 'seurat_clusters'){
        seuratObj_@meta.data$seurat_clusters <- seuratObj_@meta.data[,cluster_colname_]
    }
    
    homotypic.prop <- modelHomotypic(seuratObj_@meta.data$seurat_clusters)         
    nExp_poi <- round(0.075*nrow(seuratObj_@meta.data))  
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    seob_test_DoubletFinder <- doubletFinder_v3(seuratObj_,
                                            PCs = 1:10, pN = 0.25, 
                                            pK = pk_good, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    
    seob_test_DoubletFinder@meta.data %>% 
        rename(DF.classifications = contains('DF.classifications')) %>%
		select(-contains('pANN')) -> seob_test_DoubletFinder@meta.data
    
    ### density plotting
    ggplot(seob_test_DoubletFinder@meta.data, aes(x = nCount_RNA)) +
        scale_y_continuous(expand = c(0,0)) +
        geom_density(aes(group = DF.classifications, fill = DF.classifications), alpha = .5) +
        scale_fill_nejm() +
        theme_classic() +
        theme(text = element_text(size = 12, face = 'bold'), legend.position = 'top') -> p
    MyPlotsSave(p, filename = str_c(title_,'Density'), output = output_, width = 8, height = 8)
    
    ### col plotting
    seob_test_DoubletFinder@meta.data %>%
        group_by(seurat_clusters, DF.classifications) %>% count() %>%
        ggplot(aes(x = seurat_clusters, y = n)) +
            geom_col(aes(fill = DF.classifications), position = 'dodge', width = 0.8) +
            scale_y_continuous(expand = c(0,0)) +
            labs(y = 'cell number') +
            scale_fill_manual(values = c('#A9A9A9','#3CB371')) +
            theme_classic() +
            theme(text = element_text(size = 12, face = 'bold'), legend.position = 'top', 
                  axis.text.x = element_text(size = 12, angle = 90, vjust = .5, hjust = 1)) -> p1
    MyPlotsSave(p1, filename = str_c(title_,'Colplot'), output = output_, width = 8, height = 8)
    
    ### T-sne plotting
    seob_test_DoubletFinder@meta.data %>%
        rownames_to_column(var = 'barcodes') %>%
        left_join(as.data.frame(seob_test_DoubletFinder@reductions$umap@cell.embeddings) %>% 
                  rownames_to_column(var = 'barcodes'), by = 'barcodes') -> ggdata

    ggdata %>% as.data.frame() %>% 
        group_by(seurat_clusters) %>% 
        summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) -> label_posi

    ggplot(ggdata, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(color = seurat_clusters), size = .5) +
        geom_point(data = . %>% filter(DF.classifications == 'Doublet'), color = 'black', alpha = .3, size = .5) +
        geom_label(data = label_posi, aes(label = seurat_clusters), size = 1) + 
        guides(colour = guide_legend(override.aes = list(size = 5))) +
        ggtitle('Based on DoubletFinder') +
        scale_color_manual(values = getPalette(length(unique(seob_test_DoubletFinder@meta.data$seurat_clusters)))) +
        theme_classic() +
        theme(text = element_text(size = 12, face = 'bold')) -> p2
    MyPlotsSave(p2, filename = str_c(title_,'Umap'), output = output_, width = 8, height = 8)
    
    return(seob_test_DoubletFinder)
}
