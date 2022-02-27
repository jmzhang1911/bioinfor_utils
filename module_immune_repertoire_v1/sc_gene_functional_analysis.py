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
from myrunner import MyRunner, MyPath
from sc_utils import ScBasic

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class ScGeneFunctionalAnalysis(ScBasic):
    _ = Path(__file__).parent.absolute()
    enricher = _ / 'module_enrich/annotation_enrichment.pl'
    ppi_network_runner = _ / 'module_enrich/ppi_network.pl'
    tf_analysis_runner = _ / 'module_enrich/tf_analysis.pl'

    def __init__(self, config, id_list, statistic_file: str, all_cluster_marker_avg: str):
        """
        statistic_file: sc_basic_analysis中的所有*.statistic文件夹 -> 功能富集分析，ppi
            - 主要包括每个样本statistic
            -  cluster_diff_integrated
            -  group_diff_integrated -> 暂时将会被差分成group1_vs_group2.statistic实现并行
        all_cluster_marker_avg:差异分析中的All_cluster_Markergene_avgExp.xls -> 转录因子分析
            - sampleX.All_cluster_Markergene_avgExp.xls
            - cluster_diff_integrated.statistic/All_cluster_Markergene_avgExp.xls
        """
        super().__init__(config)
        self._id_list = id_list
        self._statistic_file = statistic_file
        # TODO:如果样本间差异分析需要实现并行，初步思路为每个差异分组创建一个group1_vs_group2.statistic
        self._statistic_file_list = self.mk_statistic_file_list()
        self._all_cluster_marker_avg = all_cluster_marker_avg.split(',')

    def mk_statistic_file_list(self):
        tmp_list = []
        for file in self._statistic_file.split(','):
            if Path(file).name.startswith('sample_diff_integrated'):
                group_statistic_dict = defaultdict(list)
                for file_ in Path(file).glob('*'):
                    if file_.name.endswith(('all_featuregene.xls', 'diff_featuregene.xls',
                                            'markergene_heatmap.pdf', 'markergene_heatmap.png')):
                        group_statistic_dict[file_.name.split('.')[0] + '.statistic'].append(file_)

                for group, file_list in group_statistic_dict.items():
                    MyPath.mkdir(group)
                    for i in file_list:
                        shutil.copy(i, group)
                    tmp_list.append(group)
            else:
                tmp_list.append(file)

        return tmp_list

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=7)
    def enrich_analysis(self):
        logging.info('doing enrich_analyse')
        cmd_list = []
        for file in self._statistic_file_list:
            cmd = '{} {} --indir {} --cfg {} --id_list {}'. \
                format(self.PERL, self.enricher, file, self._config, self._id_list)
            cmd_list.append(cmd)
        if not cmd_list:
            cmd_list = ['echo no file for enrich analysis']
        return cmd_list

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=1)
    def ppi_network(self):
        logging.info('doing ppi_network')
        cmd_list = []
        for file in self._statistic_file_list:
            cmd = '{} {} --idir {} --cfg {}'. \
                format(self.PERL, self.ppi_network_runner, file, self._config)
            cmd_list.append(cmd)
        if not cmd_list:
            cmd_list = ['echo no file for ppi analysis']
        return cmd_list

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=7)
    def tf_analysis(self):
        logging.info('doing tf_analysis')
        cmd_list = []
        for file in self._all_cluster_marker_avg:
            cmd = '{} {} --all {} --cfg {}'.format(self.PERL, self.tf_analysis_runner, file, self._config)
            cmd_list.append(cmd)
        if not cmd_list:
            cmd_list = ['echo no file for tf analysis']
        return cmd_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='running sc gene functional analysis')
    parser.add_argument('-c', '--config', type=str)
    parser.add_argument('-i', '--id_list', type=str)
    parser.add_argument('-s', '--statistic_file', type=str)
    parser.add_argument('-a', '--all_cluster_marker_avg', type=str)
    input_args = parser.parse_args()

    runner = ScGeneFunctionalAnalysis(input_args.config,
                                      input_args.id_list,
                                      input_args.statistic_file,
                                      input_args.all_cluster_marker_avg)

    executor = futures.ProcessPoolExecutor(max_workers=3)
    fs = [executor.submit(runner.enrich_analysis),
          executor.submit(runner.tf_analysis),
          executor.submit(runner.ppi_network)]

    while True:
        time.sleep(20)
        flag = True
        for f in fs:
            flag = flag and f.done()
        if flag:
            executor.shutdown()
            break
