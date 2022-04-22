.libPaths('/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/lib/R/library') 
library(optparse)

option_list <- list(
  make_option(c('-t', '--table'), type = 'character', help = 'table : gene_symbol \t value ...'), 
  make_option(c('-o', '--output'), type = 'character', help = 'plot output', default = 'GSEA_results'),
  make_option(c('-a', '--col_anno'), type = 'character', help = 'tell me colname of gene_symbol and value', 
              default = 'gene;avg_log2FC'))


opt <- parse_args(OptionParser(option_list = option_list))

library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(cowplot)
library(ggplotify)
library(foreach)
library(doParallel)

# msigdbr数据库预处理
# msigdbr(species = "Homo sapiens") %>% 
#   select(gs_cat, gs_subcat, gs_name, gene_symbol, gs_description) -> genesets
# 
# genesets %>%
#   rowwise() %>%
#   mutate(results_path = str_c(c(gs_cat, str_split(gs_subcat, ':')[[1]]), collapse= '/')) -> msigdbr_data
# 
# save(msigdbr_data, file = 'msigdbr_data.RData')
MyMkdir <- function(x){if(!dir.exists(x)){dir.create(x,recursive = T)}}

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


MyGSEAPlot <- function(results_term, raw_results, output='.'){
  set_id = results_term$ID 
  print(str_c('doing ', set_id))
  title = str_replace_all(set_id, '_', ' ') %>% str_to_lower()
  des = results_term$gs_description
  p_adjust = format(results_term$p.adjust, scientific = TRUE)
  output = file.path(output, results_term$results_path, set_id %>% str_to_lower())
  
  gseaplot(x = raw_results, geneSetID = set_id, by = 'runningScore') +
    labs(title = str_wrap(title, width = 80), 
         subtitle = str_c('p.adjust: ', p_adjust, ' (method: BH)'), 
         caption = str_wrap(des, width = 115)) +
    theme_test() -> p1
  
  MyPlotsSave(p1, filename = output, width = 8, height = 5)
  
  return(str('done -->', title))
  
}

run <- function(iterms, genelist, output='.'){
  print(str_c('doing -> ', unique(iterms$results_path)))
  
  # 创建此iterms的结果文件夹
  iterms %>%
    distinct(results_path) %>%
    pull(results_path) -> results_path;MyMkdir(file.path(output, results_path))
  
  # 获取此iterms的通路名称及对应的基因名称
  iterms %>%
    select(gs_name, gene_symbol) -> gesa_data_set
  
  # 获取通路注释及描述信息
  iterms %>%
    select(gs_name, gs_description, results_path) %>% distinct() -> doc
  
  
  # GSEA分析，对结果添加描述信息
  raw_res <- GSEA(geneList = genelist, TERM2GENE = gesa_data_set)
  raw_res@result %>% as.data.frame() %>% 
    left_join(doc, by = c('ID' = 'gs_name')) -> gsea_results
  
  # 保存当前结果
  output_filename <- str_c(str_replace_all(results_path, '/', '_'), '.csv')
  write.csv(gsea_results, 
            file = file.path(results_path, output_filename), 
            quote = F, row.names = F)
  
  # 开6个核心绘绘图
  gsea_results %>% split.data.frame(gsea_results$ID) -> all_tmp
  cl_ <- makeCluster(6)
  registerDoParallel(cl_)
  
  foreach(x = all_tmp, 
          .packages = c('tidyverse', 'enrichplot', 'ggplotify'),
          .export = c('MyMkdir', 'MyPlotsSave', 'MyGSEAPlot')
          ) %dopar%
    MyGSEAPlot(results_term = x, raw_results = raw_res, output = output)
  
  stopCluster(cl_)
}


if(T){
  ### 程序启动
  
  print('loading table and making genelist ...')
  
  gene_symbol <- str_split(opt$col_anno, pattern = ';')[[1]][1]
  value <- str_split(opt$col_anno, pattern = ';')[[1]][2]
  read.table(opt$table, sep = '\t', header = T) %>% 
    arrange(desc(value)) -> input_table
  genelist <- structure(input_table[[value]], names = input_table[[gene]])
  
  
  cl <- makeCluster(3)
  registerDoParallel(cl)
  load('./msigdbr_data.RData') # 加载整理后的基因集
  # 按照terms进行拆分数据库，每一个iterms是一个基因集的注释
  # 使用每一个iterms的基因集对genelist进行GSEA分析，三个iterms同时分析
  split.data.frame(msigdbr_data, msigdbr_data$results_path) -> tmplist
  foreach(x = tmplist,
          .packages = c('tidyverse', 'enrichplot',
                        'ggplotify', 'clusterProfiler',
                        'parallel', 'doParallel'),
          .export = c('MyMkdir', 'MyPlotsSave', 'MyGSEAPlot')) %dopar%
    run(iterms = x, genelist = genelist)
  stopCluster(cl)
}









table(genesets$gs_subcat)
genesets %>% 
  filter(gs_subcat == 'CP:KEGG') %>% 
  select(gs_name, gene_symbol) -> genesets
length(unique(genesets$gs_name))
load('./diff.exp.all_origin.RData')

diff.exp.all_origin %>% 
  filter(cluster == 'T cells, CD4+', p_val_adj >= 0) %>%
  arrange(desc(avg_log2FC)) -> genelist


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





