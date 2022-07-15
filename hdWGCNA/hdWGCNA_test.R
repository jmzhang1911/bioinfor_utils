data('mtcars')
m1 <- mtcars[, c('mpg', 'qsec')] # 性能
m2 <- mtcars[, c('cyl', 'disp', 'hp', 'drat', 'wt', 'vs', 'am')] # 车子的属性

library(ggcor)
library(ggsci)
library(cowplot) # 组合，主题
library(tidyverse)
#library(gapminder) # 数据
library(RColorBrewer)
library(ggsci)
library(ggrepel) # 标签
library(patchwork) # 拼图
library(scales) # 百分比转换
library(ggridges) # 山峦图
library(ggsignif) # p值
#library(gghalves) # 云雨图，保留箱线图，小提琴图的一半
library(ggforce) # 分面，局部放大
#library(ggpmisc) # 图中图
#library(ComplexUpset) # upset图

quickcor(m2, type = 'full') + # upper,lower,full
  geom_circle2() 
quickcor(m2, type = 'upper', cor.test = T) + 
  geom_square() + 
  geom_mark(size = 2.5, color = 'white') + # 加入pvale
  scale_fill_gradient2(low = "#00008B", high = "#8B0000", # 双渐变色
                       mid = "#708090", midpoint = 0)
quickcor(m2, type = 'upper', cor.test = T) + 
  geom_square() + 
  geom_mark(size = 2.5, color = 'white') + # 加入pvale
  scale_fill_gradient(low = "#87CEFA", high = "#4682B4") # 单渐变色
quickcor(m2, cor.test = T) + 
  geom_square(data = get_data(type = 'upper',
                              show.diag = F)) + 
  geom_mark(data = get_data(type = 'lower', # get_data()相当于继承
                            show.diag = F), # 不显示对角线
            size = 2.5, color = 'black') + # 加入pvale
  scale_fill_gradient2(low = "#00008B", high = "#8B0000", # 双渐变色
                       mid = "#708090", midpoint = 0) +
  geom_abline(slope = -1, intercept = 8,
              linetype = 'dashed')


correlate(m1, m2, cor.test = T) #ggcor中封装的cor函数
link_cor <- correlate(m1, m2, cor.test = T) %>% 
  as_cor_tbl()  # ggcor自带转换数据格式
head(link_cor)
quickcor(m2, type = 'upper', cor.test = T) + 
  geom_square() + 
  geom_mark(size = 2.5, color = 'white') + # 加入pvale
  scale_fill_gradient(low = "#87CEFA", high = "#4682B4") + # 单渐变色
  anno_link(data = link_cor,
            aes(color = p.value, size = r))
link_cor <- correlate(m1, m2, cor.test = T) %>% 
  as_cor_tbl() %>% select(.row.names, .col.names, r, p.value) %>%
  mutate(rb = cut(r, breaks = c(-1, -0.85, -0.7, 0.7, 0.85, 1),
                  labels = c('< -0.85', '-0.85 ~ -0.7', '-0.7 ~ 0.7', '0.7 ~ 0.85', '> 0.85'))) %>%
  mutate(pb = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c('< 0.01', '0.01 ~ 0.05', '> 0.05'))) 
quickcor(m2, type = 'upper', cor.test = T) + 
  geom_square() + 
  geom_mark(size = 2.5, color = 'white') + # 加入pvale
  scale_fill_gradient(low = "#87CEFA", high = "#4682B4") + # 单渐变色
  anno_link(data = link_cor,
            aes(color = pb, size = rb)) +
  scale_size_manual(values = c(0.5, 1.5, 2, 3.5)) +
  scale_color_manual(values = c("#00008B", "#8B0000", "#708090"))
