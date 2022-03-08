#!/usr/bin/env python
from myrunner import MyRunner, MyPath
from sc_utils import ScBasic
from pathlib import Path
import argparse
import logging

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class ScBasicAnalysis(ScBasic):
    _ = Path(__file__).parent.absolute() / 'ScBasicAnalysis'
    cell_filter_runner = _ / 'SingleFilterAndMergeData.R'
    basic_analyse_runner = _ / 'SingleAnalysisTsneUmap.R'
    basic_analyse_inte_runner = _ / 'samplesIntegrated.R'
    basic_analyse_inte_runner2 = _ / 'SamplesDiffInCluster.R'
    basic_analyse_inte_runner3 = _ / 'SamplesDiffInGroup.R'

    avg = 'All_cluster_Markergene_avgExp.xls'

    def __init__(self, config, sample_cellranger_result: str):
        super().__init__(config)
        self._sample_cellranger_result = sample_cellranger_result
        self._sample_file = ''
        self.get_sc_samples()

    def get_sc_samples(self):
        sc_sample_list = []
        for res in self._sample_cellranger_result.split(','):
            if res.endswith('sc.Result'):
                sc_sample_list.append(res)
        self._sample_cellranger_result = ','.join(sc_sample_list)

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=1)
    def cell_filter(self):
        logging.info('------> ScBasicAnalysis is running cell_filter')
        cmd = '{} {} -i {} -o filtered -u {} -g {} -P {} -c {}'.format(
            self.RSCRIPT,
            self.cell_filter_runner,
            self._sample_cellranger_result,  # -i
            self._config_dict['minUMI'],
            self._config_dict['minGene'],
            self._config_dict['maxpct'],
            self._config_dict['mincell'])
        cmd += ' && cp -r filtered/symbol.list ./'

        return [cmd]

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(cmd_name='sample_diff_inte_analyse', threads_num=3)
    def sample_diff_inte_analyse(self):
        """output : 1) clusterDiff_integrated 2) sampleDiff_integrated"""
        logging.info('------> ScBasicAnalysis is running sample_diff_inte_analyse')

        if len(self._sample_file) == 1:
            cmd_list = 'echo there is only one sample, skipping inte analysis'
            return [cmd_list]

        # 1) clusterDiff_integrated
        MyPath.mkdir('cluster_diff_integrated')
        sample_diff_cluster_cmd = '{} {} -R analysed_integrated/single_seruat.Rds -o cluster_diff_integrated ' \
                                  '-f {} -m {} -n {} -i {}'. \
            format(self.RSCRIPT,
                   self.basic_analyse_inte_runner2,
                   self._config_dict['fold'],
                   self._config_dict['minexp'],
                   self._config_dict['topn'],
                   self._config_dict['Symbol'])

        diff_cmd = self._extra_args(sample_diff_cluster_cmd)
        diff_cmd += ' && cp -r cluster_diff_integrated/statistic cluster_diff_integrated.statistic'
        diff_cmd += ' && cp -r cluster_diff_integrated/statistic/{} cluster_diff_integrated.{}'. \
            format(self.avg, self.avg)

        # 整合数据clusters之间的差异分析先跑，脚本脚本内占用10个线程
        MyRunner.runner(cmd_list=[diff_cmd], threads_num=1)

        # 2) sampleDiff_integrated
        MyPath.mkdir('sample_diff_integrated')
        cmd_list = []
        if self._diff_dict:
            for group, vs in self._diff_dict.items():
                if not self._config_dict.get('FDR', None):
                    logging.warning('please check config file for : FDR or pvalue')
                    break

                sample_diff_group_cmd = '{} {} -R analysed_integrated/single_seruat.Rds -o sample_diff_integrated ' \
                                        '-g {} -i {} -m {} -f {}'. \
                    format(self.RSCRIPT,
                           self.basic_analyse_inte_runner3, vs,
                           self._config_dict['Symbol'],
                           self._config_dict['minexp'],
                           self._config_dict['fold'])

                if self._config_dict.get('FDR', None):
                    sample_diff_group_cmd += ' -q {}'.format(self._config_dict['FDR'])

                if self._config_dict.get('pvalue', None):
                    sample_diff_group_cmd += ' -P {}'.format(self._config_dict['pvalue'])

                sample_diff_group_cmd += ' && cp -r sample_diff_integrated/statistic sample_diff_integrated.statistic'
                cmd_list.append(sample_diff_group_cmd)

        if not cmd_list:
            cmd_list = ['there is no diff groups, skipping ...']

        return cmd_list

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=1)
    def basic_analyse_inte(self):
        """
        output: 1) analysed_integrated
        """
        if len(self._sample_file) == 1:
            inte_cmd = 'echo there is only one sample, skipping inte analysis'
            return [inte_cmd]

        logging.info('------> ScBasicAnalysis is running basic_analyse_inte')
        MyPath.mkdir('analysed_integrated')

        inte_cmd = '{} {} -d filtered/upload.Rdata -R filtered/upload.Rdata -r {} -o {} -i {}'. \
            format(self.RSCRIPT,
                   self.basic_analyse_inte_runner,
                   self._config_dict['resolution'],
                   'analysed_integrated',
                   self._config_dict['Symbol'])

        if self._config_dict['SCT'] == 'true':
            inte_cmd += ' --SCT'

        inte_cmd += '&& cp -r analysed_integrated/{} analysed_integrated.{}'.format(self.avg, self.avg)
        inte_cmd += '&& cp -r analysed_integrated/single_seruat.Rds analysed_integrated.single_seruat.Rds'
        inte_cmd += '&& cp -r analysed_integrated/single_analysis.Rda analysed_integrated.single_analysis.Rda'

        return [inte_cmd]

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=5)
    def basic_analyse(self):
        logging.info('------> ScBasicAnalysis is running basic_analyse')
        MyPath.mkdir('analysed')
        sample_file = [rds.absolute() for rds in Path('filtered/singleSample/').rglob('*.Rds')]
        self._sample_file = sample_file

        cmd_list = []

        for sample in sample_file:
            basic_cmd = '{} {} -R {} -o {} -a {} -i {} -f {} -r {} -m {} -n {}'. \
                format(self.RSCRIPT,
                       self.basic_analyse_runner,
                       sample,  # -R
                       Path('analysed') / sample.parent.name.split('.')[0],  # -o
                       sample.parent.name.split('.')[0],  # -a from sample1_sc.sc.Result, to set sample name
                       self._config_dict['Symbol'],
                       self._config_dict['fold'],
                       self._config_dict['resolution'],
                       self._config_dict['minexp'],
                       self._config_dict['topn'])

            basic_cmd = self._extra_args(basic_cmd)
            # 将一些后续分析文件拷贝至 wd
            basic_cmd += ' && cp -r analysed/{}/clusterDiff {}.statistic'. \
                format(sample.parent.name, sample.parent.name)
            basic_cmd += ' && cp -r analysed/{}/single_seruat.Rds {}.single_seruat.Rds'. \
                format(sample.parent.name, sample.parent.name)
            basic_cmd += ' && cp -r analysed/{}/single_analysis.Rda {}.single_analysis.Rda'. \
                format(sample.parent.name, sample.parent.name)

            # TODO：警报信息
            basic_cmd += ' && cp -r analysed/{}/{}.{} ./'. \
                format(sample.parent.name, sample.parent.name, self.avg)

            cmd_list.append(basic_cmd)

        return cmd_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='running sc basic analysis')
    parser.add_argument('-c', '--config', type=str)
    parser.add_argument('-s', '--sample_cellranger_result', type=str)
    input_args = parser.parse_args()

    runner = ScBasicAnalysis(input_args.config, input_args.sample_cellranger_result)

    runner.cell_filter()
    runner.basic_analyse()
    runner.basic_analyse_inte()
    runner.sample_diff_inte_analyse()
