.libPaths('/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/lib/R/library')
library(optparse)
option_list <- list(
  make_option(c('--total'), type = 'character'),
  make_option(c('--proportions'), type = 'character'),
  make_option(c('--filename'), type = 'character'),
  make_option(c('--output'), type = 'character'))
opt <- parse_args(OptionParser(option_list = option_list))


library(tidyverse)
library(ggrepel)
library(ggsci)
library(cowplot)
library(patchwork)

proportions_plot <- function(data){
  print(str_c('processing ', basename(data)))
  read.table(data, sep = ',', header = T) %>% rename(cells = X) %>%
    mutate(cells = as.character(cells)) %>%
    pivot_longer(cols = -1, names_to = 'class', values_to = 'value') %>%
    group_by(cells) %>% arrange(class) %>%
    mutate(pos = cumsum(value),
           class = factor(class, levels = c('unspliced','spliced'))) %>%
    ggplot(aes(x = value, y = cells)) +
    geom_col(position = 'fill', aes(fill = class)) +
    geom_text(color = 'white', size = 3,
              aes(y = cells, 
                  x = pos - value * 0.5, 
                  label = scales::percent(value, accuracy = 0.2))) +
    scale_fill_jama() +
    labs(y = '', x = 'proportions') +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = 'top',
          text = element_text(size = 15, face = 'bold')) -> p
  return(p)
}

pie_plot <- function(data){
  print(str_c('processing ', basename(data)))
  read.table(data, header = T, sep = ',') %>% select(-X) %>%
    pivot_longer(cols = spliced:unspliced, 
                 values_to = 'value', names_to = 'class') %>% arrange(class) %>%
    mutate(pos = cumsum(value), 
           class = factor(class, levels = c('unspliced','spliced'))) %>%
    ggplot(aes(x = '', y = value)) +
    geom_col(aes(fill = class)) +
    coord_polar(theta = 'y') +
    geom_text(aes(y = pos - 0.5 * value, 
                  label = str_c(class, '\n', value * 100, '%')),
              color = 'white') + 
    scale_fill_jama() +
    theme_nothing() -> p
  
  return(p)
}

pie_p = pie_plot(opt$total)
proportions_p = proportions_plot(opt$proportions)
res <- pie_p | proportions_p

ggsave(plot = res, filename = str_c(opt$filename, '_proportions.png'),
       path = opt$output, width = 12, height = 6, dpi = 1200)
