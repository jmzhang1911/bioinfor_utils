library(optparse)


option_list <- list(
  make_option(c('-s', '--seob_obj'), type = 'character', help = 'Seurat Obj'), 
  make_option(c('-o', '--output'), type = 'character', help = 'output dir', default = 'hdWGCNA_results'),
  make_option(c('-c', '--config'), type = 'character', help = 'config file')
)


opt <- parse_args(OptionParser(option_list = option_list))


suppressMessages({
  library(funr)
  library(hdWGCNA)
  library(tidyverse)
  library(WGCNA)
  library(Seurat)
  library(cowplot)
  library(patchwork)
  library(ggcor)
  library(igraph)
  library(circlize)
  library(ComplexHeatmap)

  source(file.path(dirname(sys.script()), '../Utils.R'))

  # using the cowplot theme for ggplot
  theme_set(theme_cowplot())

  # set random seed for reproducibility
  set.seed(12345)

  doParallel::registerDoParallel(cl = 6)
})


read_cfg <- function(cfg_file, seurat_obj){
  read.table(cfg_file, skip = '#',  sep = '\t') %>%
    rowwise() %>%
    mutate(across(.cols = everything(), 
                  .fns =  function(x){str_trim(x)})) -> df
  
  cfg_list <- list()
  for(i in df$V1){
    value <- df[df$V1 == i,]$V2
    cfg_list[i] <- ifelse(i == 'analysis_cells', str_split(value, '@')[[1]], value)
  }
  
  if(cfg_list$analysis_cells == 'all'){
    cfg_list$analysis_cells <- as.character(unique(seurat_obj@meta.data[[cfg_list$column_cell_type]]))
  }
  
  cfg_list$assays <- ifelse('RNA' %in% names(seurat_obj@assays), 'RNA', 'Spatial')
  sample_info <- length(unique(seurat_obj@meta.data[[cfg_list$column_sample]]))
  cfg_list$sample_info <- sample_info
  
  # 写出
  assays_info <- str_c('assays', cfg_list$assays, '\n', sep = '\t')
  sample_info <- str_c('sample_number', 
                       sample_info, 
                       '\n',
                       sep = '\t')
  
  cat('\n\n#>>> temp\n', file = cfg_file, append = T)
  cat(assays_info, file = cfg_file, append = T)
  cat(sample_info, file = cfg_file, append = T)
  
  print(cfg_list)
  return(cfg_list)  
}


co_expression_network_analysis <- function(seurat_obj, cfg_list, output, cfg_file){
  # 目前hdWGCNA本身不支持Spatial，手动替换
  if(cfg_list$assays == 'Spatial'){
    names(seurat_obj@assays)[1] <- 'RNA'  
  }

  # 结果存放在 step1_co_expression_network
  output <- file.path(output, 'step1_co_expression_network');MyMkdir(output)
  res <- file.path(output, 'hdWGCNA_net_seob.rds')
  if(file.exists(res)){
    print('hdWGCNAnet exists, reading ...')
    seurat_obj <- readRDS(res)
    return(seurat_obj)
  }

  print('SetupForWGCNA ...')
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    # the gene selection approach
    gene_select = "fraction",
    # fraction of cells that a gene needs to be expressed in order to be included     
    fraction = 0.05,
    # the name of the hdWGCNA experiment             
    wgcna_name = "hdWGCNA"
  )

  print('make MetacellsByGroups ...')
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    # specify the columns in seurat_obj@meta.data to group by
    group.by = c(cfg_list$column_cell_type, cfg_list$column_sample),
    # nearest-neighbors parameter
    k = 25,
    # maximum number of shared cells between two metacells
    max_shared = 10,
    # set the Idents of the metacell seurat object
    ident.group = cfg_list$column_cell_type,
    verbose = T
  )

  # normalize metacell expression matrix:
  seurat_obj <- NormalizeMetacells(seurat_obj)

  seurat_obj <- SetDatExpr(
    seurat_obj,
    # the name of the group of interest in the group.by column
    group_name = cfg_list$analysis_cell,
    # the metadata column containing the cell type info. 
    # This same column should have also been used in MetacellsByGroups
    group.by = cfg_list$column_cell_type,
  )

  print('TestSoftPowers ...')
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    setDatExpr = FALSE, # set this to FALSE since we did this above
  )
  

  # plot the results:
  plot_list <- PlotSoftPowers(seurat_obj)

  # assemble with patchwork
  wrap_plots(plot_list, ncol=2) -> p
  MyPlotsSave(p, filename = 'TestSoftPowers', output=output, width=12, height=8)

  print('ConstructNetwork ...')
  seurat_obj <- ConstructNetwork(
    seurat_obj,
    setDatExpr=FALSE
  )

  png(file.path(output, 'Dendrogram.png'), width=800, height=400)
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
  dev.off()

  pdf(file.path(output, 'Dendrogram.pdf'), width=12, height=6)
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
  dev.off()


  # need to run ScaleData first or else harmony throws an error:
  seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

  # compute all MEs in the full single-cell dataset
  print('ModuleEigengenes ...')
  if(cfg_list$sample_info == 1){
    seurat_obj <- ModuleEigengenes(seurat_obj)
  }else{
    seurat_obj <- ModuleEigengenes(
      seurat_obj,
      group.by.vars = cfg_list$column_sample
      )
  }

  # harmonized module eigengenes:
  hMEs <- GetMEs(seurat_obj)
  write.table(head(hMEs) %>% rownames_to_column('barcode'),
              file = file.path(output, 'hMEs.xls'), 
              sep = '\t', 
              quote = F,
              row.names = F)

  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)

  # compute eigengene-based connectivity (kME):
  print('ModuleConnectivity')
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = cfg_list$column_cell_type,
    group_name = cfg_list$analysis_cell
    )
  
  hub_df_top100 <- GetHubGenes(seurat_obj, n_hubs = 100)
  hub_df_top2 <- GetHubGenes(seurat_obj, n_hubs = 2)
  write.table(hub_df_top100, file = file.path(output, 'hub_df_top100.xls'), sep = '\t', quote = F, row.names = F)
  write.table(hub_df_top2, file = file.path(output, 'hub_df_top2.xls'), sep = '\t', quote = F, row.names = F)
  
  cat(str_c('power_threshold', seurat_obj@misc$hdWGCNA$wgcna_params$power, '\n', sep = '\t'), file = cfg_file, append = T)
  print('saving hdWGCNA_net_seob ...')
  saveRDS(seurat_obj, file.path(output, 'hdWGCNA_net_seob.rds'))
  print('done network analysis')
  
  return(seurat_obj)
}


net_visualize <- function(seurat_obj, cfg_list, output){
  # 可视化模块结果存放 step2_module_plots
  print('net_visualize ...')
  TOM <- './TOM/all_ConsensusTOM-block.1.rda'
  if(!file.exists(TOM)){
    MyMkdir('TOM')
    file.copy('./ConsensusTOM-block.1.rda', TOM)
  }
  
  output <- file.path(output, 'step2_module_plots');MyMkdir(output)
  
  # get hub genes
  p <- PlotKMEs(seurat_obj, ncol=5)
  MyPlotsSave(p, filename = 'kME', output=output, width=12, height=8)

  hub_df <- GetHubGenes(seurat_obj, n_hubs = 50) # top50 hubgene
  write.table(hub_df, file=file.path(output, 'hub_df_top50.xls'), sep = '\t')

  # make a featureplot of hMEs for each module
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features='hMEs', # plot the hMEs
    order=TRUE # order so the points with highest hMEs are on top
    )
  wrap_plots(plot_list, ncol=6) -> p
  MyPlotsSave(p, filename = 'kME_umap', output=output, width=12, height=8)

  # stitch together with patchwork
  wrap_plots(plot_list, ncol=6)

  # plot module correlagram
  ModuleCorrelogram(seurat_obj) -> df
  png(file.path(output, 'Correlogram.png'), width=1000, height=1000)
  corrplot::corrplot(df$corr, type = 'upper', 
                     method = 'square',
                     tl.col = 'black',
                     tl.srt = 45,
                     tl.cex = 1.4)
  dev.off()

  pdf(file.path(output, 'Correlogram.pdf'), width=10, height=10)
  corrplot::corrplot(df$corr, type = 'upper',
                     method = 'square', 
                     tl.col = 'black', 
                     tl.srt = 45, 
                     tl.cex = 1.4)
  dev.off()

  # Volin
  MEs <- GetMEs(seurat_obj, harmonized=TRUE)
  mods <- colnames(MEs); mods <- mods[mods != 'grey']
  # add hMEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  MyMkdir(file.path(output, 'vlnplot'))
  lapply(unique(mods), function(x){
    VlnPlot(seurat_obj, features = x, 
            group.by = cfg_list$column_cell_type, pt.size = 0) + 
      geom_boxplot(width=.25, fill='white') + 
      xlab('') + ylab('hME') + NoLegend() -> p
    MyPlotsSave(p, filename = str_c(x,'_VlnPlot'), 
                output=file.path(output, 'vlnplot'), width=8, height=6)
    })

  # Individual module network plots
  ModuleNetworkPlot(seurat_obj, outdir=file.path(output, 'ModuleNetworkPlot'))

  # hubgene network
  png(file.path(output, 'HubGeneNetworkPlot.png'), width=6000, height=6000, res = 300)
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other=5,
    edge_prop = 0.75,
    mods = 'all',
    vertex.label.cex = 1.2,
    hub.vertex.size = 6,
    other.vertex.size = 1.5,
    edge.alpha = 0.4
  )
  dev.off()

  pdf(file.path(output, 'HubGeneNetworkPlot.pdf'), width=20, height=20)
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other=5,
    edge_prop = 0.75,
    mods = 'all',
    vertex.label.cex = 1.2,
    hub.vertex.size = 6,
    other.vertex.size = 1.5,
    edge.alpha = 0.4
  )
  dev.off()
  
  # HUB umap
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10, # number of hub genes to include for the UMAP embedding
    n_neighbors=15, # neighbors parameter for UMAP
    min_dist=0.1 # min distance between points in UMAP space
  )
  umap_df <- GetModuleUMAP(seurat_obj)
  umap_df %>%
    ggplot(aes(x=UMAP1, y=UMAP2)) +
      geom_point(
        color=umap_df$color, # color each point by WGCNA module
        size=umap_df$kME*2 # size of each point based on intramodular connectivity
      ) +
      umap_theme() -> p
  MyPlotsSave(p, filename = 'UMAP_hub', output=output, width=8, height=8)

  png(file.path(output, 'net_UMAP_hub.png'), width=4000, height=4000, res=300)
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1, # proportion of edges to sample (20% here)
    label_hubs=2 ,# how many hub genes to plot per module?
    keep_grey_edges=FALSE
  )
  dev.off()
  
  pdf(file.path(output, 'net_UMAP_hub.pdf'), width=10, height=10)
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1, # proportion of edges to sample (20% here)
    label_hubs=2 ,# how many hub genes to plot per module?
    keep_grey_edges=FALSE
  )
  dev.off()
 
  return(seurat_obj) 
}

module_trait_correlation <- function(seurat_obj, cfg_list, output){
  # 模块与性状的相关性 step3_module_trait_correlation
  output <- file.path(output, 'step3_module_trait_correlation');MyMkdir(output)
  
  end = ncol(seurat_obj@meta.data)
  seurat_obj@meta.data %>% rownames_to_column('barcodes') %>% 
    pivot_wider(names_from=!!as.symbol(cfg_list$column_cell_type), 
                values_from=!!as.symbol(cfg_list$column_cell_type)) %>% 
    mutate(across(.col=-c(1:end), .fns=function(x)(ifelse(is.na(x), 0, 1))),
           across(.col=-c(1:end), .fns=as.factor)) %>%
    column_to_rownames('barcodes') -> seurat_obj@meta.data
  
  
  if('nCount_RNA' %in% colnames(seurat_obj@meta.data)){
    tmp_traits <- 'nCount_RNA'
  }else{
    tmp_traits <- NA
  }
  
  if('nFeature_RNA' %in% colnames(seurat_obj@meta.data)){
    if(is.na(tmp_traits)){
      tmp_traits <- 'nFeature_RNA'
    }else{
      tmp_traits <- c(tmp_traits, 'nFeature_RNA')
    }
  }else{
    tmp_traits <- NA
  }
  
  if(is.null(cfg_list$traits)){
    if(is.na(tmp_traits)){
      tmp_traits <- cfg_list$analysis_cells
    }else{
      tmp_traits <- c(tmp_traits, cfg_list$analysis_cells)
    }
  }else{
    if(is.na(tmp_traits)){
      tmp_traits <- c(cfg_list$analysis_cells, cfg_list$traits)
    }else{
      tmp_traits <- c(tmp_traits, cfg_list$analysis_cells, cfg_list$traits)
    }
  }
  
  print(tmp_traits)
  seurat_obj <- ModuleTraitCorrelation(
    seurat_obj,
    traits = tmp_traits
  )
  
  mt_cor <- GetModuleTraitCorrelation(seurat_obj)
  cor_df <- t(mt_cor$cor$all_cells)
  fdr_df <- t(mt_cor$fdr$all_cells)
  
  write.table(cor_df, file = file.path(output, 'cor_df.xls'), sep = '\t')
  write.table(fdr_df, file = file.path(output, 'fdr_df.xls'), sep = '\t')
  
  row_df_anno <- data.frame(module = rownames(cor_df), color = rownames(cor_df))
  rownames(row_df_anno) <- rownames(cor_df)
  
  module_col_fun <- row_df_anno$module
  names(module_col_fun) <- row_df_anno$module
  
  col_fun = colorRamp2(c(-1, 0, 1), c("#87CEFA", "white", "#CC2121"))

  Heatmap(cor_df,
          border = TRUE,
          col = col_fun,
          name = 'correlation', 
          show_row_names = T,
          column_names_rot = 45, column_title = 'module trait correlation',
          right_annotation = rowAnnotation(modules =  row_df_anno[['color']],
                                           col = list(modules = module_col_fun),
                                           show_legend = FALSE),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", fdr_df[i, j]), x, y, gp = gpar(fontsize = 10))
          },
          width = 12, height = 8) -> p1
  
  Heatmap(cor_df,
          border = TRUE,
          col = col_fun,
          name = 'correlation', 
          show_row_names = T,
          column_names_rot = 45, column_title = 'module trait correlation',
          right_annotation = rowAnnotation(modules =  row_df_anno[['color']],
                                           col = list(modules = module_col_fun),
                                           show_legend = FALSE),
          cell_fun = function(j, i, x, y, width, height, fill) {
            text = fdr_df[i, j]
            if(fdr_df[i, j] < 0.001){
              text = '***'
            }else if(fdr_df[i, j] < 0.01){
              text = '**'
            }else if(fdr_df[i, j] < 0.05){
              text = '**'
            }else{
              text = ''
            }
            grid.text(text, x, y, gp = gpar(fontsize = 10))
          },
          width = 12, height = 8) -> p2
  
  pdf(file.path(output, '1_module_trait_correlation.pdf'), width = 12, height = 8)
  draw(p1)
  dev.off()
  
  png(file.path(output, '1_module_trait_correlation.png'), height = 1200, width = 1800, res = 200)
  draw(p1)
  dev.off()
  
  pdf(file.path(output, '2_module_trait_correlation.pdf'), width = 12, height = 8)
  draw(p2)
  dev.off()
  
  png(file.path(output, '2_module_trait_correlation.png'), height = 1200, width = 1800, res = 200)
  draw(p2)
  dev.off()
}

if(T){
  print('reading rds ...')
  seurat_obj <- readRDS(opt$seob_obj)
  cfg_list <- read_cfg(opt$config, seurat_obj)
  seurat_obj <- co_expression_network_analysis(seurat_obj, cfg_list, opt$output, opt$config)
  seurat_obj <- net_visualize(seurat_obj, cfg_list, opt$output)
  seurat_obj <- module_trait_correlation(seurat_obj, cfg_list, opt$output)
}



