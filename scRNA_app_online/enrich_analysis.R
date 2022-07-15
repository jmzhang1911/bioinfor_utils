
suppressMessages(library(optparse))

option_list <- list(
  make_option(c('-i', '--table'), type = 'character', help = 'input table'),
  make_option(c('-c', '--gene_colname'), type = 'character', help = 'colname of gene symbol', default = 'gene'),
  make_option(c('-g', '--GO.info'), type = 'character', help = 'GO.info'),
  make_option(c('-k', '--KEGG.info'), type = 'character', help = 'KEGG.info'),
  make_option(c('-o', '--output'), type = 'character', help = 'output', default = 'GO_KEGG_enrichment_results'),
  make_option(c('-t', '--title'), type = 'character', help = 'title', default = 'DE genes')
)
opt <- parse_args(OptionParser(option_list = option_list))

suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(ggplotify))

MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}}

SavePlot <- function(data, od, filename, width = 6, height = 4, scale = .85){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  data = as.ggplot(data, scale = scale)
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, dpi = 400, bg = 'white')
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, bg = 'white')
}

GetGeneRatio <- function(x){
  tmp <- as.numeric(str_split(x, '/')[[1]])
  return(tmp[1]/tmp[2])
}

GoEnrich <- function(genelist, info, title='DE genes', output='GO_KEGG_enrichment_results'){
  MyMkdir(output)
  TERM2GENE = select(info, V1, V3)
  TERM2NAME = select(info, V1, V2)

  tmp <- table(genelist %in% info$V3)
  print(str_c(tmp[2], ' genes in GO.info and ', tmp[1], ' not in GO.info'))
  enricher(genelist, 
           TERM2GENE = TERM2GENE, 
           TERM2NAME = TERM2NAME,
           pvalueCutoff=1,
           qvalueCutoff=1,
           pAdjustMethod = "fdr") %>% 
    as.data.frame() %>% 
    left_join(select(GO, V1, Ontology = V4) %>% distinct(), by = c('ID'='V1')) -> df
  
  if(dim(df)[1] < 1){
    return(FALSE)
  }
  
  write.table(df, file = file.path(output, 'GO_results.txt'), sep = '\t', quote = F, row.names = F)
  
  df %>%
    #filter(p.adjust < 0.1) %>% 
    group_by(Ontology) %>%
    arrange(Ontology, desc(Count)) %>% 
    slice_head(n = 8) %>%
    mutate(Description = str_wrap(Description, width = 60),
           Description = factor(Description, levels = unique(Description))) -> plot_data
  
  plot_data %>% 
    ggplot(aes(x = Description, y = -log10(pvalue))) +
    geom_col(aes(fill = Ontology), width = .8) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = '', title = str_c(title, ' GO enrichment')) + 
    scale_fill_npg() +
    theme_classic() +
    theme(axis.text.x = element_text(family="serif", 
                                     angle=315, size = 6,
                                     hjust=0, vjust=.8,
                                     face = "plain"),
          text = element_text(family="serif"),
          axis.text.y = element_text(size = 10),
          title = element_text(size = 10, face = 'bold'),
          legend.position = 'bottom',
          legend.spacing.x = unit(1.0, 'cm')) -> p
    SavePlot(data = p, od = file.path(output), filename = 'GO_col_plot', width = 13, height = 8)
  
  plot_data %>% 
    rowwise() %>%
    mutate(GeneRatio = GetGeneRatio(GeneRatio)) %>%
    ggplot(aes(x = Description, y = GeneRatio)) +
    geom_point(aes(color = pvalue, size = Count))+
    scale_size_continuous(range  = c(1, 5)) +
    facet_wrap(~Ontology, scales = 'free', nrow = 1) +
    scale_colour_gradient(low = "OrangeRed", high = "MidnightBlue") + 
    labs(x = '', title = str_c(title, ' GO enrichment')) + 
    theme_minimal() +
    theme(axis.text.x = element_text(family="serif", 
                                     angle=315, size = 8,
                                     hjust=0, vjust=.8,
                                     face = "plain"),
          text = element_text(family="serif", size = 8), 
          axis.text.y = element_text(size = 10),
          title = element_text(size = 10, face = 'bold'),
          strip.text = element_text(size = 10, face = 'bold', family="serif"),
          legend.position = 'bottom',
          legend.spacing.x = unit(1.0, 'cm')) -> p
    SavePlot(data = p, od = file.path(output), filename = 'GO_point_plot', width = 13, height = 8)
}

KeggEnrich <- function(genelist, info, title='DE genes', output='GO_KEGG_enrichment_results'){
  tmp <- table(genelist %in% info$V3)
  print(str_c(tmp[2], ' genes in KEGG.info and ', tmp[1], ' not in KEGG.info'))
  
  TERM2GENE = select(info, V1, V3)
  TERM2NAME = select(info, V1, V2)
  enricher(genelist, 
           TERM2GENE = TERM2GENE, 
           TERM2NAME = TERM2NAME, 
           pvalueCutoff=1,
           qvalueCutoff=1,
           pAdjustMethod = "fdr") %>% as.data.frame() -> df
  
  if(dim(df)[1] < 1){
    return(FALSE)
  }
  
  write.table(df, file = file.path(output, 'KEGG_results.txt'), sep = '\t', quote = F, row.names = F)
  
  df %>% 
    rowwise() %>%
    mutate(GeneRatio = GetGeneRatio(GeneRatio)) %>%
    filter(p.adjust < 0.1) %>% ungroup() %>%
    arrange(desc(GeneRatio)) %>%
    slice_head(n = 35) %>%
    mutate(Description = str_wrap(Description, width = 60),
           Description = factor(Description, levels = unique(Description))) -> plot_data
  if(nrow(plot_data) == 0){
    return()
  }
  
  plot_data %>% 
    ggplot(aes(x = Description, y = GeneRatio)) +
    geom_point(aes(color = pvalue, size = Count))+
    scale_size_continuous(range  = c(1, 4)) +
    scale_colour_gradient(low = "OrangeRed", high = "MidnightBlue") + 
    labs(x = '', title = str_c(title, ' KEGG enrichment')) + 
    theme_minimal() +
    theme(axis.text.x = element_text(family="serif", 
                                     angle=315, size = 6,
                                     hjust=0, vjust=.8,
                                     face = "plain"),
          text = element_text(family="serif", size = 10), 
          axis.text.y = element_text(size = 10),
          title = element_text(size = 10, face = 'bold'),
          legend.position = 'bottom',
          legend.spacing.x = unit(1.0, 'cm')) -> p
  
  SavePlot(data = p, od = file.path(output), filename = 'KEGG_point_plot', width = 13, height = 10)
}


marker.gene <- read.csv(opt$table, sep="\t", header=T)[,opt$gene_colname]
print('doing enrich analysis ...')
if(!is.null(opt$GO.info)){
  GO <- read.csv(opt$GO.info, header=F, sep = '\t') %>% distinct()
  print('doing GO analysis')
  GoEnrich(genelist = marker.gene, info = GO, output = opt$output, title = opt$title)
}

if(!is.null(opt$KEGG.info)){
  KEGG <- read.csv(opt$KEGG.info, header=F, sep = '\t') %>% distinct()
  print('doing KEGG analysis')
  KeggEnrich(genelist = marker.gene, info = KEGG, output = opt$output, title = opt$title)
}
