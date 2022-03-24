#!/usr/bin/env python

from myrunner import MyRunner, MyPath
from pathlib import Path
import scanpy as sc
import loompy as lp
import numpy as np
import argparse
import logging

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class PyScenic:
    _ = Path(__file__).parent
    PYTHON = '/share/nas1/zhangjm/software/miniconda3/envs/pyscenic/bin/python'
    RSCRIPT = '/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/bin/Rscript'
    RUN_pyscenic = '/share/nas1/zhangjm/software/miniconda3/envs/pyscenic/bin/pyscenic'
    pp = _ / 'pp.R'

    def __init__(self, seurat_obj, species, cell_type='cellType', output='pySCENIC_results', threads=20):
        """
        threads: 如果loom文件大小超过0.5G程序无法运行，建议对seurat进行筛选
        """
        self.output = output
        self.seurat_obj = seurat_obj
        self.cell_type = cell_type
        self.database_config = self.get_anno_data(species)
        self.threads = threads
        self.loom = ''

        if Path(self.loom).stat().st_size >= 500 * (10 ** 6):
            raise Exception('seurat obj is too big to run pySCENIC ... please filter it!!!')

        MyPath.mkdir(output)

    def get_anno_data(self, species):
        """
        按照物种增加配置文件
        - Motif_annotation_database
        - Ranking_database
        """
        motif_anno_database = self._ / 'Motif_annotation_database'
        rank_database = self._ / 'Ranking_database'
        database_config = {}
        if species == 'human':
            database_config['motif_anno'] = motif_anno_database / 'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
            database_config['rank_database'] = rank_database / 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'
            database_config['TF_list'] = self._ / 'TF_list/hs_hgnc_tfs.txt'
        else:
            raise Exception('wrong species, only support human')

        return database_config

    def run_pp(self):
        """
        传入seurat导出RNA count的loom格式，并对loom文件进行修改（解决后续BUG问题）
        """
        logging.info('getting seurat@RNA$count matrix and convert it to loom file ...')

        cmd = '{} {} --seurat_Obj {} --output {}'.format(self.RSCRIPT, self.pp, self.seurat_obj, self.output)
        MyRunner.runner([cmd], threads_num=1)

        # pyscenic后续分析会遇到bug，解决思路是修改这里的loom文件
        logging.info('reforming loom file ...')
        results = self.output + '/seob_obj.loom'
        if Path(results).exists():
            adata = sc.read_loom(Path(results))
            row_attrs = {
                "Gene": np.array(adata.var.index),
            }
            col_attrs = {
                "CellID": np.array(adata.obs.index),
                "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
                "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
            }

            for key, values in adata.obs.items():
                col_attrs[key] = np.array(values)

            loom = '{}/seob_obj_reformed.loom'.format(self.output)
            lp.create(loom, adata.X.transpose(), row_attrs, col_attrs)
            self.loom = loom

        else:
            raise Exception('there is no seob_obj_reformed.loom, wrong in pp.R')

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=1)
    def run_grn(self):
        logging.info('doing the step1 grn analysis ...')

        res = '{}/adj.sample.tsv'.format(self.output)
        if Path(res).exists():
            cmd = 'echo grn analysis has benn done, skipping ...'
        else:
            cmd = '{} {} grn {} {} --num_workers {} --output {}'. \
                format(self.PYTHON,
                       self.RUN_pyscenic,
                       self.loom,
                       self.database_config['TF_list'],
                       self.threads,
                       res)
            
        return [cmd]


if __name__ == '__main__':
    desc = """
    Version: Version v1.0
    Contact: zhangjm <zhangjm@biomarker.com.cn>
    Program Date: 2022.03.24
    Description: pySCENIC analysis
    """

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--seurat_obj', type=str, help='seurat_obj')
    parser.add_argument('-c', '--cell_type', type=str, help='colnames of cell type')
    parser.add_argument('-o', '--output', type=str, help='output dirname', default='pySCENIC_results')
    parser.add_argument('-t', '--threads', type=str, help='threads', default='20')
    parser.add_argument('-p', '--species', type=str, help='species:[human]')
    input_args = parser.parse_args()

    ps = PyScenic(seurat_obj=input_args.seurat_obj,
                  threads=input_args.threads,
                  species=input_args.species,
                  output=input_args.output,
                  cell_type=input_args.cell_type)
    ps.run_pp()
    ps.run_grn()
