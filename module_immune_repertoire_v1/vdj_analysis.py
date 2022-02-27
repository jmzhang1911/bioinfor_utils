#!/usr/bin/env python
from collections import defaultdict
from concurrent import futures
from pathlib import Path
import argparse
import logging
import shutil
import time
import sys

sys.path.append(str(Path(__file__).parent))
from myrunner import MyPath, MyRunner
from sc_utils import ScBasic

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class VdjAnalysis(ScBasic):
    _ = Path(__file__).parent.absolute()
    BCR_runner = _ / 'VDJ/BCR.R'
    TCR_runner = _ / 'VDJ/TCR.R'
    Stat_runner = _ / 'VDJ/Stat.R'
    scRepertoire_runner = Path(__file__).parent / 'VDJ/scRepertoire.R'

    def __init__(self, sample_cellranger_result: str, config, rda_files):
        """
        sample_cellranger_result:
            - sample_cellranger_result: sample1-sc.sc.Results,sample2-sc.sc.Results,sample1-t.t.Results ...
            - 需要获取其中t/b样本中中的sample1-t.t.filtered_contig_annotations.csv
        config:配置文件detail.cfg
        rda_files:单细胞分析结果中所有的rda文件
        """
        super().__init__(config)
        # TODO:需要改变cellranger的输出，及每个样本的sample.filtered_contig_annotations
        self._annotations = sample_cellranger_result.split(',')
        self._rda_files = rda_files.split(',')

        self._sample_seurat_rda = ''
        self._inte_seurat_rda = ''

        self._get_anno_dict()  # 解析所有的.csv文件
        self._split_rda()  # 解析所有的.rda文件

    def _split_rda(self):
        """对所有的rda文件进行分类：单细胞结果中inte.rda和sample.rda"""
        logging.info('splitting rds files into single_seruat and single.Rds')
        sample_rda = []
        for file in self._rda_files:
            if str(Path(file).name).startswith('analysed_integrated'):
                self._inte_seurat_rda = file
            else:
                sample_rda.append(file)
        self._sample_seurat_rda = ','.join(sample_rda)

    def _get_anno_dict(self):
        """读取所有样本的filtered_contig_annotations.csv，基于tb进行分类"""
        tmp_dict = defaultdict(list)
        for _ in [Path(file) for file in self._annotations if not str(file).endswith('sc.Result')]:
            for anno in _.glob('*.filtered_contig_annotations.csv'):
                sample, tb, _ = Path(anno).stem.split('.')
                tmp_dict[tb].append((sample, str(anno)))

        self._anno_dict = tmp_dict

    def _sc_repertoire(self, cmd, input_, output_, type_):
        input_r = self._inte_seurat_rda if Path(self._inte_seurat_rda).exists() else self._sample_seurat_rda
        cmd += '{} && {} {} -i {} -o {} -r {} -t {}'. \
            format(self.ENV, self.RSCRIPT, self.scRepertoire_runner, input_, output_, input_r, type_)
        return cmd

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper()
    def t_analysis(self):
        logging.info('doing t_analysis')
        MyPath.mkdir('t_results/input', 't_results/combined')
        for sample, file in self._anno_dict['t']:
            shutil.copy(file, 't_results/input/' + sample + '.csv')

        cmd = '{} {} -i t_results/input/ -o {} -t TRA -s {}'. \
            format(self.RSCRIPT, self.TCR_runner,
                   't_results/TRA',
                   self._config_dict['Species'])

        cmd += ' && {} {} -i t_results/input/ -o {} -t TRB -s {}'. \
            format(self.RSCRIPT, self.TCR_runner,
                   't_results/TRB',
                   self._config_dict['Species'])

        cmd += ' && {} {} -i t_results/input/ -o t_results/stat -t TCR'. \
            format(self.RSCRIPT, self.Stat_runner)

        cmd = self._sc_repertoire(cmd, 't_results/input', 't_results/combined', 'T')

        if not cmd:
            return ['echo no T file, skipping t analysis']

        return [cmd]

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper()
    def b_analysis(self):
        logging.info('doing b_analysis')
        MyPath.mkdir('b_results/input', 'b_results/combined')
        for sample, file in self._anno_dict['b']:
            shutil.copy(file, 'b_results/input/' + sample + '.csv')

        cmd = '{} {} -i b_results/input/ -o {} -t IGH -s {}'. \
            format(self.RSCRIPT, self.BCR_runner,
                   'b_results/IGH',
                   self._config_dict['Species'])

        cmd += ' && {} {} -i b_results/input/ -o {} -t IGK -s {}'. \
            format(self.RSCRIPT, self.BCR_runner,
                   'b_results/IGK',
                   self._config_dict['Species'])

        cmd += ' && {} {} -i b_results/input/ -o {} -t IGL -s {}'. \
            format(self.RSCRIPT, self.BCR_runner,
                   'b_results/IGL',
                   self._config_dict['Species'])

        cmd += ' && {} {} -i b_results/input/ -o b_results/stat -t BCR'. \
            format(self.RSCRIPT, self.Stat_runner)

        cmd = self._sc_repertoire(cmd, 'b_results/input', 'b_results/combined', 'B')

        if not cmd:
            return ['echo no T file, skipping t analysis']
        
        return [cmd]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='running vdj basic analysis')
    parser.add_argument('-f', '--sample_cellranger_result', type=str)
    parser.add_argument('-c', '--config', type=str)
    parser.add_argument('-i', '--rda_files', type=str)
    input_args = parser.parse_args()
    runner = VdjAnalysis(input_args.sample_cellranger_result,
                         input_args.config,
                         input_args.rda_files)

    executor = futures.ProcessPoolExecutor(max_workers=2)
    fs = [executor.submit(runner.t_analysis), executor.submit(runner.b_analysis)]

    while True:
        time.sleep(20)
        flag = True
        for f in fs:
            flag = flag and f.done()
        if flag:
            executor.shutdown()
            break
