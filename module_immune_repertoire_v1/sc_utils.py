#!/usr/bin/env python
from myrunner import MyRunner
from pathlib import Path


class ScBasic:
    # Python
    python2 = 'python2'
    python3 = 'python3'

    # Rscript
    RSCRIPT = 'Rscript'

    # Perl
    PERL = 'perl'

    # Env：R依赖的动态库
    ENV = ' && export LD_LIBRARY_PATH={}/2.2.1/lib/:$LD_LIBRARY_PATH && ' \
          'export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/iconv/1.16/lib:$LD_LIBRARY_PATH && ' \
          'export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/Anaconda3/2019.03/lib/:$LD_LIBRARY_PATH && ' \
          'export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/gcc/5.5.0/lib:/share/nas2/genome/biosoft/gcc/5.5.0/lib64:$LD_LIBRARY_PATH && ' \
          'export PATH=/share/nas2/genome/biosoft/gcc/5.5.0/bin/:$PATH'.format(Path(__file__).parent / 'gsl')

    def __init__(self, config, data_config=''):
        self._config = config  # detail.config文件
        self._data_config = data_config  # data.config文件

        # detail.config dict 文件
        self._config_dict = self._get_config()

        # 差异分组dict，{group1: sample1_sample2_vs_sampleA_sampleB, group2:a_b_vs_1_2}
        self._diff_dict = {k: samples for k, samples in self._config_dict.items() if 'group' in k and '-sc' in samples}

        # 所有的执行命令存入list
        self._cmd_list = []

    def _get_config(self):
        with open(self._config, 'r') as f:
            tmp_dict = {}
            for line in f:
                if len(line.split()) == 1:
                    tmp_dict[line.strip()] = ''
                else:
                    k, v = line.split()
                    tmp_dict[k] = v.strip()

            return tmp_dict

    def _extra_args(self, basic_cmd):
        """一些差异分析常用的额外参数"""
        if self._config_dict['SCT'] == 'true':
            basic_cmd += ' --SCT'

        if self._config_dict.get('FDR', None):
            basic_cmd += ' -q {}'.format(self._config_dict['FDR'])

        if self._config_dict.get('pvalue', None):
            basic_cmd += ' -p {}'.format(self._config_dict['pvalue'])

        return basic_cmd

    @MyRunner.count_running_time
    def run_all_cmd(self):
        MyRunner.runner(self._cmd_list)
        self._cmd_list = []
        return self._cmd_list
