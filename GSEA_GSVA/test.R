library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(cowplot)
library(Seurat)
library(ggplotify)
library(foreach)
library(doParallel)


MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}else(cat('existed dir\n'))}

MyPlotsSave <- function(plot, filename='__tmp', output='./', width=10, height=10, dpi = 500){
  # Give me a plot, save it as ggplot object
  ggsave(as.ggplot(plot, scale = 0.95), 
         filename = str_c(filename, '.pdf'), 
         path = output, width = width, 
         height = height, 
         limitsize = FALSE)
  
  width = if_else(width>35, 20, width)
  height = if_else(width>35, 20, height)
  ggsave(as.ggplot(plot, scale = 0.95), 
         filename = str_c(filename, '.png'),
         path = output, 
         width = width, 
         height = height, 
         dpi = dpi)
}


# msigdbr数据库预处理
# msigdbr(species = "Homo sapiens") %>% 
#   select(gs_cat, gs_subcat, gs_name, gene_symbol, gs_description) -> genesets
# 
# genesets %>%
#   rowwise() %>%
#   mutate(results_path = str_c(c(gs_cat, str_split(gs_subcat, ':')[[1]]), collapse= '/')) -> msigdbr_data
# 
# save(msigdbr_data, file = 'msigdbr_data.RData')



load('./msigdbr_data.RData')

MyGSEAPlot <- function(results_term, raw_results){
  set_id = results_term$ID 
  print(str_c('doing ', set_id))
  title = str_replace_all(set_id, '_', ' ') %>% str_to_lower()
  des = results_term$gs_description
  p_adjust = format(results_term$p.adjust, scientific = TRUE)
  output = file.path(results_term$results_path, set_id %>% str_to_lower())
  
  gseaplot(x = raw_results, geneSetID = set_id, by = 'runningScore') +
    labs(title = str_wrap(title, width = 80), 
         subtitle = str_c('p.adjust: ', p_adjust, ' (method: BH)'), 
         caption = str_wrap(des, width = 115)) +
    theme_test() -> p1
  
  MyPlotsSave(p1, filename = output, width = 8, height = 5)
  
  return(str('done -->', title))
  
}

run <- function(iterms, genelist){
  print(str_c('doing -> ', unique(iterms$results_path)))
  
  iterms %>%
    distinct(results_path) %>%
    pull(results_path) -> results_path;MyMkdir(results_path)
  
  iterms %>%
    select(gs_name, gene_symbol) -> gesa_data_set
  
  iterms %>%
    select(gs_name, gs_description, results_path) %>% distinct() -> doc
  
  
  raw_res <- GSEA(geneList = genelist, TERM2GENE = gesa_data_set)
  raw_res@result %>% as.data.frame() %>% 
    left_join(doc, by = c('ID' = 'gs_name')) -> gsea_results
  
  # 保存当前结果
  output_filename <- str_c(str_replace_all(results_path, '/', '_'), '.csv')
  write.csv(gsea_results, 
            file = file.path(results_path, output_filename), 
            quote = F, row.names = F)
  
  gsea_results %>% split.data.frame(gsea_results$ID) -> all_tmp
  cl_ <- makeCluster(6)
  registerDoParallel(cl_)
  
  foreach(x = all_tmp, 
          .packages = c('tidyverse', 'enrichplot', 'ggplotify'),
          .export = c('MyMkdir', 'MyPlotsSave', 'MyGSEAPlot')
          ) %dopar%
    MyGSEAPlot(results_term = x, raw_results = raw_res)
  
  stopCluster(cl_)
}


cl <- makeCluster(3)
registerDoParallel(cl)
split.data.frame(msigdbr_data, msigdbr_data$results_path) -> tmplist
foreach(x = tmplist,
        .packages = c('tidyverse', 'enrichplot',
                      'ggplotify', 'clusterProfiler',
                      'parallel', 'doParallel'),
        .export = c('MyMkdir', 'MyPlotsSave', 'MyGSEAPlot')) %dopar%
  run(iterms = x, genelist = genelist)
stopCluster(cl)





table(genesets$gs_subcat)
genesets %>% 
  filter(gs_subcat == 'CP:KEGG') %>% 
  select(gs_name, gene_symbol) -> genesets
length(unique(genesets$gs_name))
load('./diff.exp.all_origin.RData')

diff.exp.all_origin %>% 
  filter(cluster == 'T cells, CD4+', p_val_adj >= 0) %>%
  arrange(desc(avg_log2FC)) -> genelist
genelist <- structure(genelist$avg_log2FC, names = genelist$gene)

res <- GSEA(geneList = genelist, TERM2GENE = genesets)
gseaplot(res, geneSetID = 'KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY', 
         title = 'KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY', by = 'runningScore') -> p

ca <- 'Genes having at least one occurence of the motif AAACCAC in their 3 untranslated region. The motif represents putative target (that is, seed match) of human mature miRNA hsa-miR-140 (v7.1 miRBase).'
p +
  labs(title = str_wrap('KEGG T CELL RECEPTOR SIGNALING PATHWAY', width = 40), 
       subtitle = 'p.adjust: 0.05 (method: BH)', 
       caption = str_wrap(ca, width = 130)) +
  theme_test() -> p1
  
ggsave(p1, filename = 'tmp.png', width = 8, height = 6, dpi = 500)

a <- genesets$gene_symbol
split(genesets, a) %>% 
  lapply(function(x){return(x$gene_symbol)}) -> gene_set_gsva





