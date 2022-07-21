#!/usr/bin/env python
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent.parent))

from bioinfor_tools.cmd_runner import CmdRunner
from bioinfor_tools.utils import BioUtils
import argparse
import logging


class Cpdb(BioUtils):
    PYTHON = '/share/nas1/zhangjm/software/miniconda3/envs/cpdb/bin/python3'
    RSCRIPT = '/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/bin/Rscript'
    cpdb = '/share/nas1/zhangjm/software/miniconda3/envs/cpdb/bin/cellphonedb'

    def __init__(self, seurat_obj, species='human', result_dir='cpdb_results', cell_type_col='cellType', threads=15,
                 group_config='None'):

        super().__init__(module=__file__)
        self.cellphonedb_pre = self._module_dir / 'cellphonedb_pre.R'

        self.seurat_obj = seurat_obj
        self.species = species
        self.cell_type_col = cell_type_col
        self.result_dir = result_dir
        self.threads = threads
        self.group_config = group_config

    @CmdRunner.cmd_wrapper(use_qsub=True)
    def run_get_matrix(self):
        """获取meta_data以及gene_count"""

        self.mkdir(self.result_dir)
        logging.info('running get matrix')
        cmd = '{} {} --MyMakeMatrix --seurat_Obj {} --species {} --cell_type {} --results {} --group_config {}'. \
            format(self.RSCRIPT, self.cellphonedb_pre,
                   self.seurat_obj, self.species,
                   self.cell_type_col, self.result_dir,
                   self.group_config)

        return [cmd]

    @CmdRunner.cmd_wrapper(use_qsub=True)
    def run_cellphone_db(self):
        """运行cellphonedb，计算，点图，及热图"""
        logging.info('running cellphonedb')
        file_list = list(Path(self.result_dir).glob('*_cpdb'))
        cmd_list = []
        for file in file_list:
            # 计算
            cmd = '{} {} method statistical_analysis {} {} --threads {} --output-path {}'. \
                format(self.PYTHON,
                       self.cpdb,
                       file.absolute() / 'meta_data.txt',
                       file.absolute() / 'gene_count.txt', self.threads,
                       Path(self.result_dir) / file.name / 'out')
            # 点图
            cmd += '&& {} {} plot dot_plot --means-path {} --pvalues-path {} --output-path {}'. \
                format(self.PYTHON,
                       self.cpdb,
                       Path(self.result_dir) / file.name / 'out/means.txt',
                       Path(self.result_dir) / file.name / 'out/pvalues.txt',
                       Path(self.result_dir) / file.name / 'out')

            # 热图，并生成net_count.txt用于可视化
            cmd += ' && {} {} plot heatmap_plot  {} --pvalues-path {} --output-path {}'. \
                format(self.PYTHON,
                       self.cpdb,
                       file.absolute() / 'meta_data.txt',
                       Path(self.result_dir) / file.name / 'out/pvalues.txt',
                       Path(self.result_dir) / file.name / 'out')

            # 将几个pdf转为png
            cmd += ' && convert {0}/plot.pdf {0}/plot.png' \
                   ' && convert {0}/heatmap_count.pdf {0}/heatmap_count.png' \
                   ' && convert {0}/heatmap_log_count.pdf {0}/heatmap_log_count.png'. \
                format(*(str(Path(self.result_dir) / file.name / 'out'),))

            cmd_list.append(cmd)

        return cmd_list

    @CmdRunner.cmd_wrapper(use_qsub=True)
    def run_plot(self):
        """使用igraph绘制网络图"""
        logging.info('running plotting')
        file_list = list(Path(self.result_dir).glob('*_cpdb'))
        cmd_list = []
        for file in file_list:
            cmd = '{} {} --MyCpdbPlot --count_network {}/count_network.txt --output {}'. \
                format(self.RSCRIPT, self.cellphonedb_pre,
                       Path(self.result_dir) / file.name / 'out',
                       Path(self.result_dir) / file.name / 'out')
            cmd_list.append(cmd)

        return cmd_list

    @CmdRunner.cmd_wrapper(use_qsub=True)
    def run_output(self):
        """打包结果文件，添加readme文件"""
        self.copy_readme(self.result_dir)
        cmd = 'tar -zcvf {}.tar.gz {} --exclude gene_count.txt --exclude meta_data.txt'. \
            format(self.result_dir,
                   Path(self.result_dir).name,
                   Path(self.result_dir).name)

        return [cmd]


if __name__ == '__main__':
    doc_path_group_config = Path(__file__).parent / 'group_config.cfg'

    desc = """
    Version: Version v1.1
    Contact: zhangjm <zhangjm@biomarker.com.cn>
    Program Date: 2021.12.01
    Program UpDate: 2022.02.18
    Description: cellphonedb analysis
     - Support human, mouse, rat
     - Support custom setting groups conveniently
     - using data instead of RNA
    """

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--seurat_obj', type=str, help='Seuratobj')
    parser.add_argument('-p', '--species', type=str, default='human', help='species [human,mouse,rat]')
    parser.add_argument('-o', '--result_dir', type=str, default='cpdb_results', help='output dir')
    parser.add_argument('-c', '--cell_type_col', type=str, default='cellType', help='colname of celltype')
    parser.add_argument('-t', '--threads', type=int, default=15, help='number of threads')
    parser.add_argument('-g', '--group_config', type=str, default='None',
                        help='config file ps:{}, if not set will do all the seurat_obj whole together'.format(
                            doc_path_group_config))
    input_args = parser.parse_args()

    runner = Cpdb(seurat_obj=input_args.seurat_obj, species=input_args.species,
                  result_dir=input_args.result_dir, cell_type_col=input_args.cell_type_col,
                  threads=input_args.threads, group_config=input_args.group_config)

    runner.make_summary()
    runner.run_get_matrix()
    runner.run_cellphone_db()
    runner.run_plot()
    runner.run_output()
    runner.make_summary()
