library(msigdbr)
library(tidyverse)


msigdbr(species = "Mus musculus") %>%
  select(gs_cat, gs_subcat, gs_name, gene_symbol, gs_description) -> genesets

genesets %>%
  rowwise() %>%
  mutate(results_path = str_c(c(gs_cat, str_split(gs_subcat, ':')[[1]]), collapse= '/')) -> msigdbr_data

save(msigdbr_data, file = 'mm_msigdbr_data.RData')
