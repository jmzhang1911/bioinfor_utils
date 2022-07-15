#!/usr/bin/env python
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))
from myrunner import MyRunner, MyPath, make_summary
import argparse
import logging
import shutil

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class MyInferCnv:
    RSCRIPT = '/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/bin/Rscript'
    script_path = Path(__file__).parent
    infercnv_pre_processing_runner = script_path / 'infercnv_pre_processing.R'
    run_infercnv_runner = script_path / 'run_infercnv.R'
    infercnv_subclusters_runner = script_path / 'infercnv_subclusters.R'
    read_me = script_path / 'infercvn_readme.zip'

    debug_cmd = 'ulimit -s unlimited'

    def __init__(self, seob, group_config, cell_config, modes, gene_order, thresholds='auto',
                 output='InferCNV_results'):
        self.seob = seob
        self.group_config = group_config
        self.cell_config = cell_config
        self.modes = modes
        self.gene_order = gene_order
        self.thresholds = thresholds
        self.output = output

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper()
    def pp(self):
        if Path('step1_infercnv_pp/myinfer_obj.RData').exists():
            cmd = 'echo step1_infercnv_pp done!'
            return [cmd]

        cmd = '{} {} --seurat_Obj {} --group_config {} --cell_config {}'. \
            format(self.RSCRIPT,
                   self.infercnv_pre_processing_runner,
                   self.seob,
                   self.group_config,
                   self.cell_config)
        return [cmd]

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper()
    def run_infercnv(self):
        if Path('step2_infercnv_temp').exists():
            cmd = 'echo step2_infercnv_temp done!'
            return [cmd]

        cmd = '{} && {} {} --myinfercnv_obj step1_infercnv_pp/myinfer_obj.RData --modes {} --gene_order {}'. \
            format(self.debug_cmd,
                   self.RSCRIPT,
                   self.run_infercnv_runner,
                   self.modes,
                   self.gene_order)
        return [cmd]

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper()
    def infercnv_subclusters(self):
        if Path('step3_infercnv_subclusters').exists():
            cmd = 'echo step3_infercnv_subclusters done!'
            return [cmd]

        cmd = '{} {} --myinfercnv_obj step2_infercnv_temp/myinfer_obj.RData --thresholds {}'. \
            format(self.RSCRIPT,
                   self.infercnv_subclusters_runner,
                   self.thresholds)
        return [cmd]

    @MyRunner.cmd_wrapper()
    def make_results(self):
        backup = 'BackUp'
        MyPath.mkdir(self.output, backup)

        # 备份完整结果文件
        for _ in ['step1_infercnv_pp', 'step2_infercnv_temp', 'step3_infercnv_subclusters']:
            shutil.copy(_, backup)

        # 展示文件
        step2_file_list = ['infercnv.observations.txt ', 'infercnv.references.txt',
                           'HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat ',
                           'HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat']
        shutil.copy('step1_infercnv_pp/anno_file.txt', self.output)
        for _ in step2_file_list:
            shutil.copy('step2_infercnv_temp/{}'.format(_), self.output)

        for _ in Path('step3_infercnv_subclusters').glob('**/*'):
            shutil.copy(_, self.output)
        shutil.copy(self.read_me, self.output)

        cmd_list = ['tar -zcvf BackUp.tar.gz Backup',
                    'tar -zcvf {}.tar.gz {}'.format(self.output, self.output)]

        return cmd_list


if __name__ == '__main__':
    doc_path_cell_config = Path(__file__).parent / 'cell_config.cfg'
    doc_path_group_config = Path(__file__).parent / 'group_config.cfg'
    doc_path_gene_order_human = Path(__file__).parent / 'Homo_gene_order_file.txt'
    doc_path_gene_order_mouse = Path(__file__).parent / 'Mus_musculus_GRCm38_release95_gene_order_file.txt'

    desc = """
    Version: Version v1.0
    Contact: zhangjm <zhangjm@biomarker.com.cn>
    Program Date: 2022.02.18
    Description: infercnv analysis
     - Support subcluster mode and cell mode
     - Support custom setting groups, setting reference and observing cell groups conveniently
    """

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--seurat_obj', type=str, help='Seuratobj')
    parser.add_argument('-g', '--group_config', type=str,
                        help='group_config ps:{}'.format(doc_path_group_config.absolute()))
    parser.add_argument('-c', '--cell_config', type=str,
                        help='cell_config ps:{}'.format(doc_path_cell_config.absolute()))
    parser.add_argument('-m', '--modes', type=str,
                        help='modes [subclusters|cells] default: subclusters', default='subclusters')
    parser.add_argument('-r', '--gene_order', type=str,
                        help='gene_order reference:\nhuman:{}\nmouse:{}'.format(doc_path_gene_order_human,
                                                                                doc_path_gene_order_mouse))
    parser.add_argument('-t', '--thresholds', type=str,
                        help='thresholds for color if auto is not perfect, ps 0.8,1,1.2 [default auto]',
                        default='auto')
    input_args = parser.parse_args()

    infer_cnv = MyInferCnv(seob=input_args.seurat_obj,
                           group_config=input_args.group_config,
                           cell_config=input_args.cell_config,
                           modes=input_args.modes,
                           gene_order=input_args.gene_order,
                           thresholds=input_args.thresholds)

    make_summary(Path(__file__), status='doing')
    infer_cnv.pp()
    infer_cnv.run_infercnv()
    infer_cnv.infercnv_subclusters()
    make_summary(Path(__file__), status='done')
