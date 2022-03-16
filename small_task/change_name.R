
library(tidyverse)

gene_id2trans_id <- read.table('./gene_id2trans_id.txt', sep = '\t') %>%
  distinct() %>% 
  mutate(V2 = str_replace_all(V2, ':', '_'))
symbol_list <- read.table('./symbol.list')

MyChangeKegg <- function(x){
  all_gene <- str_split(x, ';') %>% unlist() %>%
    str_remove_all('.gene')
  
  res <- data.frame(origin = all_gene) %>% 
    left_join(gene_id2trans_id, by = c('origin'='V2')) %>% drop_na() %>% pull(V1) 
  
  if(length(res)>1){
    res %<>%str_c(collapse = ';') %>% str_c(';')
  }
  
  return(res)
}

Known_longest_transcript_fa_Kegg <- read_delim("Known.longest_transcript.fa.Kegg.pathway", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE) 
Known_longest_transcript_fa_Kegg %>%
  rowwise() %>% 
  mutate(Gene_id = as.character(Gene_id),
         Gene_id = MyChangeKegg(x = Gene_id)) -> kegg_pathway

Known_longest_transcript_fa_GO_list <- read_delim("Known.longest_transcript.fa.GO.list.txt", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  col_names = FALSE, trim_ws = TRUE)
Known_longest_transcript_fa_GO_list %>%
  rowwise() %>% 
  mutate(X1 = MyChangeKegg(x = X1)) -> go_list


Known_longest_transcript_fa_GO_anno <- read_delim("Known.longest_transcript.fa.GO.anno.txt", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  trim_ws = TRUE)

Known_longest_transcript_fa_GO_anno %>% 
  rowwise() %>%
  mutate(`#Gene` = MyChangeKegg(x = `#Gene`)) -> go_anno


ppi %>% mutate(`#Query_id1` = str_remove(`#Query_id1`, '.gene'), 
               `Query_id2` = str_remove(`Query_id2`, '.gene')) %>% 
  left_join(gene_id2trans_id, by = c('#Query_id1'='V2')) %>% 
  mutate(`#Query_id1` = V1) %>% select(-V1) %>%
  left_join(gene_id2trans_id, by = c('Query_id2'='V2')) %>%
  mutate(Query_id2 = V1) %>% select(-V1) -> res

# save
dir.create('results')


kegg_pathway %>%
  write.table(file = file.path('results', 'Known.longest_transcript.fa.Kegg.pathway'), 
              sep = '\t', row.names = F, col.names = T, quote = F)

go_list %>%
  write.table(file = file.path('results', 'Known.longest_transcript.fa.GO.list.txt'), 
              sep = '\t', row.names = F, col.names = T, quote = F)

go_anno %>%
  write.table(file = file.path('results', 'Known.longest_transcript.fa.GO.anno.txt'), 
              sep = '\t', row.names = F, col.names = T, quote = F)



