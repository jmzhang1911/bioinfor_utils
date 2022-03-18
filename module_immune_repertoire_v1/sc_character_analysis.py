#!/usr/bin/env python
from myrunner import MyPath, MyRunner
from concurrent import futures
from sc_utils import ScBasic
from pathlib import Path
import argparse
import logging
import time

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class ScCharacterAnalysis(ScBasic):
    _ = Path(__file__).parent.absolute()
    trace_analyse_runner = _ / 'ScBasicAnalysis/TraceAnalysis.R'
    cell_anno_analyse_runner = _ / 'ScBasicAnalysis/cellAnnotation.R'
    cell_cycle_analyse_runner = _ / 'cell_cycle/cell_cycle.R'

    def __init__(self, config, rds: str):
        """
        rds: sample1.single_seruat.Rds,sample2.single_seruat.Rds,analysed_integrated.single_seruat.Rds'
        """
        super().__init__(config)
        self._rds = rds  # 所有的rds文件
        self._rds_dict = self._make_rds_dict()

    def _make_rds_dict(self):
        logging.info('splitting rds files into single_seruat and single.Rds')
        tmp_dict = {}
        for file in self._rds.split(','):
            sample = Path(file).name.replace('.single_seruat.Rds', '')
            tmp_dict[sample] = file
        return tmp_dict

    def _extra_args(self, basic_cmd):
        """
        一些差异分析常用的额外参数
        与sc_utils._extra_args不一致，需要单独再写一遍了，真无语！！！
        在此处FDR对应为f
        """
        if self._config_dict['SCT'] == 'true':
            basic_cmd += ' --SCT'

        if self._config_dict.get('FDR', None):
            basic_cmd += ' -f {}'.format(self._config_dict['FDR'])

        if self._config_dict.get('pvalue', None):
            basic_cmd += ' -p {}'.format(self._config_dict['pvalue'])

        return basic_cmd

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=2)
    def trace_analyse(self):
        logging.info('doing the trace_analyse')
        cmd_list = []
        for sample, rds in self._rds_dict.items():
            cmd = '{} {} -S {} -o {}.cell_trace -a {} -e {} -C {}'. \
                format(self.RSCRIPT, self.trace_analyse_runner, rds, sample, sample,
                       self._config_dict['minexp'],
                       self._config_dict['mincell'])
            cmd = self._extra_args(cmd)
            cmd_list.append(cmd)
        return cmd_list

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=2)
    def cell_anno_analyse(self):
        logging.info('doing the cell_anno_analyse')
        if self._config_dict['refData'] == 'null':
            return ['echo no refData has benn choose, skipping cell anno analysis']

        cmd_list = []
        refdata = (Path(__file__).parent / 'ref.singleR' / 'ref.{}.Rds'.format(self._config_dict['refData']))
        for sample, rds in self._rds_dict.items():
            cmd = '{} {} -r {} -f {} -o {}.cell_typeAnno -s {}'. \
                format(self.RSCRIPT, self.cell_anno_analyse_runner, rds,
                       refdata, sample, sample)
            cmd_list.append(cmd)
        return cmd_list

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=2)
    def cell_cycle_analyse(self):
        logging.info('cell_cycle_analyse')
        cmd_list = []
        for sample, rds in self._rds_dict.items():
            cmd = '{} {} -s {} -m {} -o {}.cell_cycle'. \
                format(self.RSCRIPT, self.cell_cycle_analyse_runner,
                       self._config_dict['Species'], rds,
                       sample)
            cmd_list.append(cmd)
        return cmd_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='running sc character analysis')
    parser.add_argument('-c', '--config', type=str)
    parser.add_argument('-r', '--rds', type=str)
    input_args = parser.parse_args()

    runner = ScCharacterAnalysis(input_args.config, input_args.rds)

    executor = futures.ProcessPoolExecutor(max_workers=2)
    fs = [executor.submit(runner.cell_anno_analyse),
          executor.submit(runner.cell_cycle_analyse)]

    while True:
        time.sleep(20)
        flag = True
        for f in fs:
            flag = flag and f.done()
        if flag:
            executor.shutdown()
            break

    # trace资源消耗过大，放在最后
    runner.trace_analyse()
