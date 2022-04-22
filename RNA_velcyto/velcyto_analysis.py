#!/usr/bin/env python
import shutil
import sys

sys.path.append('/share/nas1/zhangjm/workspace/MyUtils')

from concurrent.futures import ThreadPoolExecutor
from make_web_report.web_report import WebReport
from myrunner import MyRunner, MyPath
from anndata import AnnData
from pathlib import Path
import seaborn as sns
import scvelo as scv
import pandas as pd
import scanpy as sc
import numpy as np
import argparse
import logging
import loompy

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


class ScveloAnalysis:
    """
    notes:
        - ipython path:/home/zhangjm/software/python3.9.6/bin/ipython
        - please using UMAP
`   ------------------------------------------------------------------------------------------
    01): running samtools and velcyto. It takes almost 12 h per sample.
    -> step01: running samtools:
    ScveloAnalysis.run_samtools('path/to/cellranger_results/SC')

    -> step02: running velcyto
    reference:
        - seurat_GTF = '/share/nas1/zhangjm/workspace/project/AH787/RNA_velocyto/genes.gtf'
        - RMSK = '/share/nas1/zhangjm/workspace/project/AH787/RNA_velocyto/Human_rmsk.gtf'
    MyScveloAnalysis.run_velcyto(cellranger_results,
                                RMSK='path/to/repeater's.gtf',
                                seurat_GTF='path/to/cellranger.gtf',
                                overwrite= False)
    ------------------------------------------------------------------------------------------
    02): RNA velcyto analysis. It takes almost 1 h per sample.
    -> step01: group_info=False if group_info.txt doesn't exist
    vel = ScveloAnalysis('path/to/cellranger/SC','group_info.txt')

    -> step02: merged_path='path/to/merged_path' if merged loom files existed
    vel.read()

    -> step03: filter data and normal analysis
    vel.filter_velcyto()

    -> step04 optional: using seurat data(contains row_names) if needed
    vel.using_seurat_meta_umap('seurat_meta_df', 'seurat_umap')

    -> step05: vel.read_dy_model('path/to/dy_results') if `run_dy_model` has been executed
    vel.run_dy_model()

    -> step06: RNA velcyto analysis
    vel.rna_velcyto()
    -----------------------------------------------------------------------------------------
    03): making Report
    vel.run_report()
    """

    WORKERS = 8
    SAMTOOLS_THREADS = 10
    DPI = 1200

    RSCRIPT = '/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/bin/Rscript'
    scvelo_plot = Path(__file__).parent / 'scvelo_plot.R'

    def __init__(self, cellranger_results: str = False, group_info=False):
        self.cluster = 'leiden'
        self.color = 'leiden'
        self._group_info = group_info
        self._loom_dict = {}

        if cellranger_results:
            self._cellranger_results = Path(cellranger_results)
            self._loom_dict = self.__get_loom_file_dict()

        if group_info:
            self._group_info_dict = self.__get_groups()
            self.merge_loom()

    @staticmethod
    @MyRunner.cmd_wrapper(threads_num=WORKERS, test=False)
    def run_samtools(cellranger_results, thread=SAMTOOLS_THREADS, overwrite=False) -> list:
        samtools = '/share/nas1/zhangjm/software/miniconda3/bin/samtools'

        def get_cmd(sample):
            # 解决bam名称问题
            input_bam = Path(sample) / 'outs/possorted_genome_bam.bam'
            # input_bam = [bam for bam in list((Path(sample) / 'outs').glob('*possorted_genome_bam.bam'))
            #              if not str(bam).startswith('cellsorted')]
            # if not len(input_bam) == 1:
            #     raise Exception('wrong bam file in sample:{}'.format(Path(sample).name))
            # else:
            #     input_bam = input_bam[0]

            # 输入bam
            output_bam = sample / 'outs/cellsorted_possorted_genome_bam.bam'
            cmd = '{} sort -@ {} -t CB -O BAM -o {} {}'.format(samtools, thread, output_bam, input_bam)
            if Path(output_bam).exists():
                cmd = 'echo {} bam file already sorted!'.format(sample) if not overwrite else cmd

            return cmd

        cmd_list = [get_cmd(i.absolute()) for i in Path(cellranger_results).iterdir()]

        return cmd_list

    @staticmethod
    @MyRunner.cmd_wrapper(threads_num=WORKERS, test=False)
    def run_velcyto(cellranger_results, RMSK, seurat_GTF, overwrite=False) -> list:
        python3 = '/home/zhangjm/software/python3.9.6/bin/python3'
        RNA_velocyto = '/home/zhangjm/software/python3.9.6/bin/velocyto'

        def get_cmd(sample):
            cmd = '{} {} run10x -m {} {}/{} {}'.format(python3, RNA_velocyto, RMSK, cellranger_results, sample.name,
                                                       seurat_GTF)
            if (sample / 'outs/velocyto' / '{}.loom'.format(sample.name)).exists():
                cmd = 'echo {} loom file already exist!'.format(sample) if not overwrite else cmd

            return cmd

        cmd_list = [get_cmd(i.absolute()) for i in Path(cellranger_results).iterdir()]

        return cmd_list

    # 合并或者与seurat同步的loom文件存放在__loom_dict中
    @property
    def loom_dict(self):
        return self._loom_dict

    # 原始的loom文件
    def __get_loom_file_dict(self):
        tmp_dict = {sample.name: sample / 'velocyto' / (sample.name + '.loom') for sample in
                    self._cellranger_results.iterdir()}
        return tmp_dict

    # 读取需要进行分组的分组信息表
    def __get_groups(self):
        with open(self._group_info, 'r') as f:
            group_dict = {}
            for i in f:
                group = i.split()[0]
                sample = i.split()[1].split(';')
                group_dict[group] = ['{}/{}/velocyto/{}.loom'.format(self._cellranger_results, i, i) for i in sample]
        return group_dict

    # 按照分组信息表进行合并merge，存放位置默认为Merged_loom
    def merge_loom(self, output='Merged_loom'):
        def combine(loom_list: list, output_file):
            logging.info('doing -> {}'.format(output_file))
            loompy.combine(loom_list, output_file=output_file)

        MyPath.mkdir(output)
        logging.info('Merging the loom files')

        with ThreadPoolExecutor(max_workers=len(self._group_info_dict)) as executor_cmd:
            for group, sample_list in self._group_info_dict.items():
                output_loom = Path(output, group + '_combined.loom')
                executor_cmd.submit(combine, sample_list, output_loom)

        logging.info('loom files are merged!!!')
        self._loom_dict = {sample: output + '/' + sample + '_combined.loom' for sample in self._group_info_dict.keys()}

    # 基于合并或者未合并的loom_file_dict读取loom文件，并输出一个loom文件的格式用于调整seurat对象的行名
    def read_loom(self, merged_path: str = False):
        logging.info('loading...')

        if merged_path:
            # 直接读入合并的loom
            logging.info('reading the merged loom files directly...')
            loom_dict = {loom.name.replace('_combined.loom', ''): sc.read(loom, sparse=True, cache=True) for loom in
                         Path(merged_path).iterdir()}
        else:
            logging.info('reading the loom files...')
            loom_dict = {sample: sc.read(loom, sparse=True, cache=True) for sample, loom in self._loom_dict.items()}

        list(loom_dict.values())[0].obs.head(10).to_csv('scanpy_meta.data.txt')
        self._loom_dict = loom_dict

    # 将seurat的降维结果替换掉scanpy中的，替换后的loom文件存放在__loom_dict中，即最终用于分析的loom文件
    def using_seurat_meta_umap(self, seurat_meta_df, seurat_umap):
        # TODO:发现规律：
        #  Lep01-sc:AAATGCCTCACATGCAx = AAACCTGAGAGTACAT-1_1
        #  如果一直是这种规律的话可以直接替换了，最终形态是直接输入一个seurat对象

        self.color = 'seurat_clusters'
        self.cluster = 'seurat_clusters'

        seurat_meta_data_raw = pd.read_csv(seurat_meta_df)
        seurat_umap_raw = pd.read_csv(seurat_umap)

        def change_meta_umap(loom_file):
            anndata_index = loom_file.obs.index.values
            meta_data = loom_file.obs.reset_index(). \
                merge(seurat_meta_data_raw, left_on="CellID", right_on="barcodes", how='left'). \
                set_index('CellID').reindex(anndata_index)

            try:
                # 解决cluster带.0的问题
                meta_data['seurat_clusters'] = meta_data['seurat_clusters'].apply(
                    lambda x: "%d" % int(x) if ~np.isnan(x) else x).astype('category')
            except TypeError:
                meta_data['seurat_clusters'] = meta_data['seurat_clusters'].astype('category')

            loom_file.obs = meta_data
            loom_file = loom_file[loom_file.obs.dropna().index, :]

            loom_file.obsm['umap'] = seurat_umap_raw.set_index('barcodes'). \
                reindex(loom_file.obs.index).to_numpy().astype('float64')

            return loom_file

        self._loom_dict = {sample: change_meta_umap(loom) for sample, loom in self._loom_dict.items()}

        logging.info('using seurat data... done!')

    def change_colors(self):
        # 可能是版本问题，出图无法使用全局设置的颜色，这里手动设置，使用Set1的颜色类型
        for sample, anndata in self._loom_dict.items():
            categories_num = len(anndata.obs.value_counts(self.cluster))
            anndata.__dict__['_uns'][self.cluster + '_colors'] = sns.color_palette("Set1", categories_num)
            self._loom_dict[sample] = anndata

    def filter_velcyto(self, results='scvelo_results'):
        def _filter(anndata):
            # 如果不适用seurat中的数据，则基于scanpy进行聚类
            scv.pp.filter_genes(anndata, min_shared_counts=20)
            scv.pp.normalize_per_cell(anndata)
            scv.pp.filter_genes_dispersion(anndata, n_top_genes=2000)
            scv.pp.log1p(anndata)
            scv.pp.filter_and_normalize(anndata, min_shared_counts=20, n_top_genes=2000)
            scv.pp.moments(anndata, n_pcs=30, n_neighbors=30)

            # 如果不使用Seurat对象中的数据则使用scanpy做降维聚类
            if self.cluster != 'seurat_clusters':
                scv.tl.umap(anndata)
                sc.tl.leiden(anndata)

            return anndata

        for sample, raw_anndata in self._loom_dict.items():
            out_put = Path(results) / sample
            MyPath.mkdir(out_put)
            self._loom_dict[sample] = _filter(raw_anndata)

    def read_dy_model(self, dy_results='dy_model_data', used_seurat=False):
        if used_seurat:
            self.color = 'seurat_clusters'
            self.cluster = 'seurat_clusters'

        if Path(dy_results).is_file():
            self._loom_dict = {Path(dy_results).name.replace('_dy_model.h5ad', ''): scv.read(dy_results)}
        else:
            self._loom_dict = {sample.name.replace('_dy_model.h5ad', ''): scv.read(sample) for sample in
                               Path(dy_results).glob('**/*') if str(sample).endswith('h5ad')}

    def run_dy_model(self, dy_results='dy_model_data', n_jobs=30):
        MyPath.mkdir(dy_results)

        def _run_dy_model(anndata: sc.AnnData):
            scv.tl.recover_dynamics(anndata, n_jobs=n_jobs)
            scv.tl.velocity(anndata, mode='dynamical')
            scv.tl.velocity_graph(anndata, n_jobs=n_jobs)
            anndata.write(Path(dy_results) / (sample + '_dy_model.h5ad'), compression='gzip')

            return anndata

        self.change_colors()
        for sample, raw_anndata in self._loom_dict.items():
            self._loom_dict[sample] = _run_dy_model(raw_anndata)

    def _save_fig(self, p, save_as):
        p.figure.savefig(str(save_as) + '.png', bbox_inches='tight', dpi=self.DPI)
        p.figure.savefig(str(save_as) + '.svg')
        logging.info('{} has been saved'.format(Path(save_as).name))

    @staticmethod
    def run_report(config):
        """生成报告, 所有东西一定得拷贝到Web_Report"""
        des = 'Web_Report/BMK_1_data'
        MyPath.mkdir(des)
        template_file = Path(__file__).parent / 'template.txt'
        for i in Path('scvelo_results').glob('**/*'):
            if str(i).endswith('.png') or str(i).endswith('.csv'):
                shutil.copy(str(i), des)

        wr = WebReport(config=config, template=template_file, report_name='RNA-velcyto_report')
        wr.config_template()
        wr.make_report()

    @staticmethod
    def get_proportions_data(andata, groupby="clusters", use_raw=True, layers=None):
        """原图proportions展示结果奇怪，这里提出数据重新绘图，此函数返回pd对象，此部分参考源码"""
        from scvelo.core import sum

        if layers is None:
            layers = ["spliced", "unspliced", "ambigious"]

        layers_keys = [key for key in layers if key in andata.layers.keys()]
        counts_layers = [sum(andata.layers[key], axis=1) for key in layers_keys]

        if use_raw:
            ikey, obs = "initial_size_", andata.obs
            counts_layers = [
                obs[ikey + layer_key] if ikey + layer_key in obs.keys() else c
                for layer_key, c in zip(layers_keys, counts_layers)
            ]

        counts_total = np.sum(counts_layers, 0)
        counts_total += counts_total == 0
        counts_layers = np.array([counts / counts_total for counts in counts_layers])

        counts_groups = dict()
        for cluster in andata.obs[groupby].cat.categories:
            counts_groups[cluster] = np.mean(
                counts_layers[:, andata.obs[groupby] == cluster], axis=1
            )

        labels = list(counts_groups.keys())
        data = np.array(list(counts_groups.values()))

        # 比例分布图数据
        df = pd.DataFrame(data).set_index(pd.Series(labels))
        df.columns = ['spliced', 'unspliced']

        # 饼图数据
        total = pd.DataFrame({'spliced': round(np.mean(counts_layers, axis=1)[0], 3),
                              'unspliced': round(np.mean(counts_layers, axis=1)[1], 3)},
                             index=['1'])

        return df, total

    # 分析的核心内容，后续如果需要添加新的分析功能请在此处添加
    def rna_velcyto(self, results='scvelo_results'):
        def _rna_velcyto(anndata: AnnData):

            # proportions plot
            p = scv.pl.proportions(anndata, groupby=self.cluster, show=False, figsize=(15, 4))
            df, total = self.get_proportions_data(anndata, groupby=self.cluster)
            df_results = out_put / (sample + '_proportions.csv')
            total_results = out_put / (sample + '_total.csv')

            df.to_csv(df_results)
            total.to_csv(total_results)
            cmd = '{} {} --total {} --proportions {} --filename {} --output {}'. \
                format(self.RSCRIPT,
                       self.scvelo_plot,
                       total_results,
                       df_results,
                       sample,
                       out_put)
            MyRunner.runner([cmd])

            if isinstance(p, list):
                # 此图片大小不规则，因此使用默认的方式保存而非self._save_fig()
                p[0].figure.savefig(out_put / (sample + '_proportions_origin.png'))
                p[0].figure.savefig(out_put / (sample + '_proportions_origin.pdf'))
            else:
                p.figure.savefig(out_put / (sample + '_proportions_origin.png'))
                p.figure.savefig(out_put / (sample + '_proportions_origin.pdf'))

            # plotting umap
            p = scv.pl.umap(anndata, color=self.cluster, show=False)
            self._save_fig(p, out_put / (sample + '_umap'))

            # embedding_stream analysis
            p = scv.pl.velocity_embedding_stream(anndata, basis='umap', color=self.cluster, show=False)
            self._save_fig(p, out_put / (sample + '_velocity_embedding_stream_dy'))

            # velocity_embedding analysis
            p = scv.pl.velocity_embedding(anndata, arrow_length=3, arrow_size=2, show=False, color=self.cluster)
            self._save_fig(p, out_put / (sample + '_velocity_embedding_dy'))

            # PAGA analysis
            scv.tl.paga(anndata, groups=self.cluster)
            p = scv.pl.paga(anndata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, show=False)
            self._save_fig(p, out_put / (sample + '_paga_dy'))

            # velocity_pseudotime analysis
            scv.tl.velocity_pseudotime(anndata)
            p = scv.pl.scatter(anndata, color='velocity_pseudotime', cmap='gnuplot', show=False)
            self._save_fig(p, out_put / (sample + '_velocity_pseudotime_dy'))

            # latent_time analysis
            scv.tl.latent_time(anndata)
            p = scv.pl.scatter(anndata, color='latent_time', color_map='gnuplot', size=80, show=False)
            self._save_fig(p, out_put / (sample + '_latent_time_dy'))

            # drive gene
            # heatmap
            top_genes = anndata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
            scv.pl.heatmap(anndata, var_names=top_genes, sortby='latent_time', n_convolve=100,
                           col_color=self.cluster, save=str(out_put / (sample + '_drive_heatmap.png')))

            anndata.var['fit_likelihood'].sort_values(ascending=False).head(300). \
                to_csv(out_put / (sample + '_top300_fit_likehood.csv'))

            # scatter
            top_genes = anndata.var['fit_likelihood'].sort_values(ascending=False).index
            scv.pl.scatter(anndata, basis=top_genes[:15], ncols=5, frameon=False,
                           color=self.cluster, save=str(out_put / (sample + '_drive_top15.png')),
                           dpi=self.DPI)

            anndata.var['fit_likelihood'].sort_values(ascending=False).head(15). \
                to_csv(out_put / (sample + '_top15_fit_likelihood.csv'))

            return anndata

        self.change_colors()
        for sample, raw_anndata in self._loom_dict.items():
            out_put = Path(results) / sample
            MyPath.mkdir(out_put)
            self._loom_dict[sample] = _rna_velcyto(raw_anndata)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RNA Velcyto analysis', formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers()

    # run_velcyto
    parser_run_velcyto = subparsers.add_parser('run_velcyto', help='run_velcyto help')
    parser_run_velcyto.add_argument('-i', '--input', type=str, help='super cellranger outputs dir')
    parser_run_velcyto.add_argument('-t', '--max_workers', dest='max_workers', type=int, help='set thread number',
                                    default=6)
    parser_run_velcyto.add_argument('--overwrite', action='store_true', help='output folder')
    parser_run_velcyto.add_argument('-r', '--RMSK', type=str, help='repeatmasker GTF')
    parser_run_velcyto.add_argument('-s', '--seurat_GTF', type=str, help='seurat GTF')
    parser_run_velcyto.add_argument('-@', '--thread', type=str, help='samtools WORKERS', default=8)

    # run_scVelo
    parser_run_scVelo = subparsers.add_parser('run_scVelo', help='run_scVelo help')
    parser_run_scVelo.add_argument('-i', '--input', type=str, help='cellranger_results')
    parser_run_scVelo.add_argument('-n', '--n_jobs', type=int, help='set thread number', default=20)
    parser_run_scVelo.add_argument('-g', '--grouped', type=str, help='grouped file', default=False)
    parser_run_scVelo.add_argument('-S', '--Seurat_data', action='store_true', help='using seurat data')
    parser_run_scVelo.add_argument('-m', '--Seurat_meta_data', type=str, help='seurat meta.data')
    parser_run_scVelo.add_argument('-u', '--Seurat_umap', type=str, help='seurat umap data')
    parser_run_scVelo.add_argument('-c', '--config', type=str, help='report config')

    # RNA velcyto analysis
    parser_rna_velcyto = subparsers.add_parser('rna_velcyto', help='add all args showed in run_velcyto and run_scVelo')
    parser_rna_velcyto.add_argument('-i', '--input', type=str, help='super cellranger outputs dir')
    parser_rna_velcyto.add_argument('-t', '--max_workers', dest='max_workers', type=int, help='set thread number',
                                    default=8)
    parser_rna_velcyto.add_argument('--overwrite', action='store_true', help='output folder')
    parser_rna_velcyto.add_argument('-r', '--RMSK', type=str, help='repeatmasker GTF')
    parser_rna_velcyto.add_argument('-s', '--seurat_GTF', type=str, help='seurat GTF')
    parser_rna_velcyto.add_argument('-@', '--thread', type=str, help='samtools WORKERS', default=8)
    parser_rna_velcyto.add_argument('-n', '--n_jobs', type=int, help='set thread number', default=20)
    parser_rna_velcyto.add_argument('-g', '--grouped', type=str, help='grouped file', default=False)
    parser_rna_velcyto.add_argument('-S', '--Seurat_data', action='store_true', help='using seurat data')
    parser_rna_velcyto.add_argument('-m', '--Seurat_meta_data', type=str, help='seurat meta.data')
    parser_rna_velcyto.add_argument('-u', '--Seurat_umap', type=str, help='seurat umap data')
    parser_rna_velcyto.add_argument('-c', '--config', type=str, help='report config')

    parser_run_velcyto.set_defaults(action=('run_velcyto', lambda: None))
    parser_run_scVelo.set_defaults(action=('run_scVelo', lambda: None))
    parser_rna_velcyto.set_defaults(action=('rna_velcyto', lambda: None))

    args = parser.parse_args()
    arg, _ = args.action


    def up_stream():
        MyScveloAnalysis.run_samtools(cellranger_results=args.input,
                                      thread=args.thread,
                                      overwrite=args.overwrite)
        MyScveloAnalysis.run_velcyto(cellranger_results=args.input,
                                     RMSK=args.RMSK,
                                     seurat_GTF=args.seurat_GTF,
                                     overwrite=args.overwrite)


    def down_steam():
        vel = MyScveloAnalysis(args.input, group_info=args.grouped)
        vel.read_loom()

        if args.Seurat_data:
            vel.using_seurat_meta_umap(seurat_meta_df=args.Seurat_meta_data,
                                       seurat_umap=args.Seurat_umap)

        vel.filter_velcyto()
        vel.run_dy_model()
        vel.rna_velcyto()

        if args.config:
            vel.run_report(config=args.config)


    class MyScveloAnalysis(ScveloAnalysis):
        try:
            WORKERS = args.max_workers
            SAMTOOLS_THREADS = args.thread
        except AttributeError:
            logging.info('ready to go...')


    if arg == 'run_velcyto':
        up_stream()

    if arg == 'run_scVelo':
        down_steam()

    if arg == 'rna_velcyto':
        up_stream()
        down_steam()
