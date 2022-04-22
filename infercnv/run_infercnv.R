library(optparse)

option_list <- list(
  make_option(c('-s', '--myinfercnv_obj'), type = 'character', help = 'myinfercnv_obj'),
  make_option(c('-m','--modes'), type = 'character', help = 'modes[subclusters|cells]', default = 'subclusters'),
  make_option(c('-g','--gene_order'), type = 'character', help = 'gene_order')
)
opt <- parse_args(OptionParser(option_list = option_list))

library(tidyverse)
library(infercnv)
options(expressions=500000)
load(opt$myinfercnv_obj)

print('creating infercnv obj')
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= myinfercnv_obj$gene_matrix,
                                    annotations_file= myinfercnv_obj$cell_anno, 
                                    gene_order_file=  opt$gene_order,
                                    ref_group_names= myinfercnv_obj$ref_cell)

print('running infercnv')
infercnv_obj_res = infercnv::run(infercnv_obj,
                                 # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 cutoff = 0.1,
                                 out_dir = 'step2_infercnv_temp',
                                 cluster_by_groups = TRUE,
                                 num_threads = 30,
                                 # 分析到亚克隆，速度较慢
                                 analysis_mode = opt$modes,
                                 denoise = TRUE,
                                 HMM = TRUE)

print('saving results')
subclusters_cell_groupings = 'step2_infercnv_temp/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings'

# 添加信息至myinfercnv_obj
# 添加infer_obj，infercnv最终结果，用于热图值
myinfercnv_obj['infer_obj'] = 'step2_infercnv_temp/run.final.infercnv_obj'
# 添加亚克隆信息
myinfercnv_obj['subclusters_cell_groupings'] = subclusters_cell_groupings
# 添加热图阈值信息
myinfercnv_obj['thresholds'] = 'step2_infercnv_temp/infercnv.21_denoised.heatmap_thresholds.txt'
# 添加cnv得分信息表
myinfercnv_obj['infercnv.observations'] = 'step2_infercnv_temp/infercnv.observations.txt'
# 添加分析模式信息
myinfercnv_obj['mode'] = opt$modes

save(myinfercnv_obj, file = 'step2_infercnv_temp/myinfer_obj.RData')



