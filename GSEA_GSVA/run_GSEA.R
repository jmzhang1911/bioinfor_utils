.libPaths('/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/lib/R/library') 
library(optparse)

option_list <- list(
  make_option(c('-t', '--table'), type = 'character', help = 'gene_symbol value ...'), 
  make_option(c('-o', '--output'), type = 'character', help = 'plot output', default = 'GSEA_results'),
  make_option(c('-a', '--col_anno'), type = 'character', 
              help = 'tell me colname of gene_symbol and value ps:gene,avg_log2FC', 
              default = 'gene,avg_log2FC'),
  make_option(c('-s', '--species'), type = 'character', help = 'species [hs|mm]', default = 'hs'),
  make_option(c('-g', '--gs_cat'), type = 'character', help = '[C1|C2|C3|C4|C5|C6|C7|C8|H], sep by ,', default = 'all')
  )


opt <- parse_args(OptionParser(option_list = option_list))

suppressMessages({
  library(msigdbr)
  library(GSVA)
  library(tidyverse)
  library(clusterProfiler)
  library(patchwork)
  library(cowplot)
  library(ggplotify)
  library(foreach)
  library(doParallel)
  library(funr)
  source(file.path(dirname(sys.script()), '../Utils.R'))
})

MyGSEAPlot <- function(results_term, raw_results, output='.'){
  set_id = results_term$ID 
  print(str_c('doing ', set_id))
  title = str_replace_all(set_id, '_', ' ') %>% str_to_lower()
  des = results_term$gs_description
  p_adjust = format(results_term$p.adjust, scientific = TRUE)
  
  NES_score = round(results_term$NES, digits = 4)
  #output = file.path(output, results_term$results_path, set_id %>% str_to_lower())
  
  gseaplot(x = raw_results, geneSetID = set_id, by = 'runningScore') +
    labs(title = str_wrap(title, width = 80), 
         subtitle = str_c('p.adjust: ', p_adjust, ' (method: BH)\nNES:', NES_score), 
         caption = str_wrap(des, width = 115)) +
    theme_test() -> p1
  
  file_name <- set_id %>% str_to_lower() %>% basename()
  MyPlotsSave(p1,
              output = file.path(output, results_term$results_path),
              filename = file_name,
              width = 8, height = 5)
  
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
  raw_res <- tryCatch({
    GSEA(geneList = genelist, TERM2GENE = gesa_data_set)
  }, error = function(e) {
    
  })
  if(is.null(raw_res)){
    return()
  }
  
  nrow_number <- nrow(as.data.frame(raw_res@result))
  if( nrow_number< 1){
    print('no results')
    return()
  }
  
  raw_res@result %>% as.data.frame() %>% 
    left_join(doc, by = c('ID' = 'gs_name')) -> gsea_results
  
  #增加一个点图
  gsea_results %>% 
    mutate(direction = if_else(NES > 0, 'up', 'down')) %>%
    group_by(direction) %>%
    slice_max(abs(NES), n = 10) %>%
    arrange(desc(NES)) %>%
    mutate(Description = str_replace_all(Description, '_', ' '),
           Description = str_wrap(Description, width = 50), 
           Description = factor(Description)) %>% 
    ggplot(aes(x = NES, y = Description)) +
    geom_col(aes(fill = p.adjust)) +
    labs(title = str_c('TOP10 NES iterms in ', basename(results_path))) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 6)) -> p
  print(output)
  print(file.path(output, results_path))
  MyPlotsSave(plot = p, 
              output = file.path(output, results_path),
              filename = str_c(str_replace_all(results_path, '/', '_'), 
                               '_col_plot'),
              height = 8, width = 8)
  
  
  # 保存当前结果
  output_filename <- str_c(str_replace_all(results_path, '/', '_'), '.csv')
  write.csv(gsea_results, 
            file = file.path(output, results_path, output_filename), 
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


make_summary(sys.script(), 'doing')
if(T){
  ### 程序启动
  
  print('loading table and making genelist ...')
  
  gene_symbol <- str_split(opt$col_anno, pattern = ',')[[1]][1]
  value <- str_split(opt$col_anno, pattern = ',')[[1]][2]
  
  print(str_c('gene_symbol colname => ', gene_symbol))
  print(str_c('value colname => ', value))
  
  read.table(opt$table, sep = '\t', header = T) -> input_table
  genelist <- sort(structure(input_table[[value]], names = input_table[[gene_symbol]]), decreasing = T)
  
  if(opt$species == 'hs'){
    load(file.path(dirname(sys.script()), 'hs_msigdbr_data.RData'))
  }else if(opt$species == 'mm'){
    load(file.path(dirname(sys.script()), 'mm_msigdbr_data.RData'))
  }else{
    stop('wrong species, only for hs|mm') 
  }
  
  # 按照terms进行拆分数据库，每一个iterms是一个基因集的注释
  # 使用每一个iterms的基因集对genelist进行GSEA分析，三个iterms同时分析
  
  df <- msigdbr_data %>% ungroup()
  if (opt$gs_cat != 'all') {
    gs_cat_vector <- str_split(opt$gs_cat, ',')[[1]]
    df <- df %>% filter(gs_cat %in% gs_cat_vector)
    print(str_c('using ', gs_cat_vector))
  }
  
  split.data.frame(df, df$results_path) -> tmplist
  
  print('doing GSEA ...')
  cl <- makeCluster(3)
  registerDoParallel(cl)
  foreach(x = tmplist,
          .packages = c('tidyverse', 'enrichplot',
                        'ggplotify', 'clusterProfiler',
                        'parallel', 'doParallel'),
          .export = c('MyMkdir', 'MyPlotsSave', 'MyGSEAPlot')) %dopar%
    run(iterms = x, genelist = genelist, output = opt$output)
  stopCluster(cl)
}

make_summary(sys.script(), 'done')