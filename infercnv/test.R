library(tidyverse)
library(ggsci)
df <- read.table('/Volumes/My_Passport/data_of_server/infercnv.observations.txt')
df2 <- read.table('/Volumes/My_Passport/data_of_server/infercnv.observation_groupings.txt')
df3 <- read.table('/Volumes/My_Passport/data_of_server/anno_file.txt', sep ='\t')
infer <- readRDS('/Volumes/My_Passport/data_of_server/run.final.infercnv_obj')

load('/Volumes/My_Passport/data_of_server/expr.RData')

df %>% rownames_to_column('genes') %>%
  pivot_longer(cols = -1, names_to = 'barcodes', values_to = 'expr') %>%
  mutate(barcodes = str_replace_all(barcodes, '\\.', '-')) %>%
  left_join(df3, by = c('barcodes'='V1')) %>% View()



MyNormalize <- function(x){
  return((x - min(x)) / (max(x) - min(x)))
}

a <- c(1,2,3,4,5,5,6)
kmeans_res <- kmeans(expr, centers = 10)

kmeans_res <- kmeans(t(as.matrix(df %>% drop_na())), centers = 6)
results <- as.data.frame(kmeans_res$cluster) %>% rename(cluster = `kmeans_res$cluster`) %>%
  rownames_to_column('barcodes') %>%
  mutate(barcodes = str_replace_all(barcodes, '\\.', '-'))


df[1:10, 1:10] %>%
  mutate(across(.cols = everything(), .fns = MyNormalize)) %>% View()

apply(df, 2, MyNormalize) %>% View()
hist(df$TCTACATGTGTCGATT.1_9)

df %>% rownames_to_column('genes') %>%
  pivot_longer(cols = -1, names_to = 'barcodes', values_to = 'expr')  %>%
  mutate(barcodes = str_replace_all(barcodes, '\\.', '-')) %>%
  group_by(barcodes) %>%
  summarise(total = sum(expr ^ 2)) %>% 
  left_join(results, by = 'barcodes') %>%
  mutate(cluster = str_c(as.character(cluster), '_cluster')) %>%
  #group_by(cluster) %>%
  #summarise(total = sum(expr)) %>%
  ggplot(aes(x = cluster, y = total)) +
  geom_violin(aes(fill = cluster)) +
  labs(y = 'CNV scores', x = '') +
  scale_fill_nejm() + 
  theme_bw() +
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = 'top') -> p1
ggsave(p1, filename = 'tmp.png', width = 12)




if(len <= 5){ 
  n.col <- len
}else if(len <= 10){
  ifelse(len %% 2 ==0, n.col <- len/2, n.col <- len/3)
}else{
  n.col <- ceiling(len/6)
}




