library(tidyverse)
library(cowplot)

df <- read.table('/Volumes/My_Passport/data_of_server/logging_summary.txt')

df %>%
  rename(date = V1, time = V2, pwd = V3, script = V4, status = V5) %>%
  rowwise() %>%
  mutate(user = str_split(pwd, '/')[[1]][4],
         #script = str_remove_all(script, 'script_id:|run_|.py|.R'),
         script = case_when( script == 'script_id:run_cell_trace.py' ~  '轨迹分析' ,
                    script == 'script_id:velcyto_analysis.py' ~ 'RNA速度',
                    script == 'script_id:run_infercnv.py' ~ 'inferCNV',
                    script == 'script_id:run_scenic.py' ~ 'pySCENIC',
                    script == 'script_id:run_GSEA.R' ~ 'GSEA',
                    script == 'script_id:run_cellphonedb.py' ~ '细胞通讯')) %>%
  filter(status == 'status:done', user != 'zhangjm') %>% 
  count(script) %>%
  ggplot(aes(x = script, y = n)) +
  scale_y_continuous(limits = c(0, 13)) +
  geom_col(fill = '#008080') +
  geom_text(aes(label=n, y=n+0.05), position=position_dodge(0.9), vjust=0) +
  labs(x = '', y = '运行次数', title = '六月份生产部门执行成功的个性化分析') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p

source('Utils.R')

MyPlotsSave(p, filename = '202206.png', height = 6, width = 6)

font_add("NotoSansSCBlack", "/Volumes/My_Passport/Blackhole/Noto_Sans_SC/NotoSansSC-Thin.otf")
font_families()
