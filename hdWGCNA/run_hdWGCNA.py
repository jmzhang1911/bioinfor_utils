import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))

from bioinfor_tools.utils import make_summary, read_cfg, MyPath
from scRNA_app_online.enrich_analysis import EnrichAnalysis
from make_web_report.web_report import WebReport
from bioinfor_tools.cmd_runner import CmdRunner
import argparse
import logging
import shutil

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class hdWGCNA:
    PYTHON = '/share/nas1/zhangjm/software/miniconda3/envs/cpdb/bin/python3'
    RSCRIPT = '/share/nas1/zhangjm/software/miniconda3/envs/R4.1.3/bin/Rscript'
    run_hdWGCNA_runner = Path(__file__).parent / 'run_hdWGCNA.R'
    template_file = Path(__file__).parent / 'template.txt'

    def __init__(self, seob_obj, output, config):
        self.seob_obj = seob_obj
        self.output = output
        self.config = config
        self.cfg_dict = ''

    @CmdRunner.cmd_wrapper(use_qsub=True)
    def run_hdWGCNA(self):
        cmd = '{} {} --seob_obj {} --output {}  --config {}'. \
            format(self.RSCRIPT, self.run_hdWGCNA_runner, self.seob_obj, self.output, self.config)
        return [cmd]

    def run_enrich_analysis(self):
        self.cfg_dict = read_cfg(self.config)
        enrich = EnrichAnalysis(diff_results_table=Path(self.output) / 'step1_co_expression_network/hub_df_top100.xls',
                                output=Path(self.output) / 'step4_enrichment',
                                GO_info=self.cfg_dict['GO_info'],
                                KEGG_info=self.cfg_dict['KEGG_info'],
                                split_by='module',
                                symbol_list=self.cfg_dict['symbol_list'],
                                gene_column='gene_name')
        enrich.run_enrichment()

    def run_report(self):
        bmk_1 = 'Web_Report/BMK_1'
        bmk_2 = 'Web_Report/BMK_2'
        bmk_3 = 'Web_Report/BMK_3'

        MyPath.mkdir(*[bmk_1, bmk_2, bmk_3])

        for file in Path(self.output).glob('**/*'):
            if file.name in ['TestSoftPowers.png', 'Dendrogram.png']:
                shutil.copy(file, bmk_1)

            if file.name in ['Correlogram.png', 'kME.png', 'kME_umap.png', 'net_UMAP_hub.png',
                             'HubGeneNetworkPlot.png', 'hMEs.xls', 'hub_df_top2.xls']:
                shutil.copy(file, bmk_2)

            if file.name.endswith('_module_trait_correlation.png'):
                shutil.copy(file, bmk_2)

            if file.name in ['GO_col_plot.png', 'GO_point_plot.png', 'KEGG_point_plot.png']:
                shutil.copy(file, '{}/{}.{}'.format(bmk_3, file.parent.name, file.name))

        MyPath.mkdir(bmk_2 + '/ModuleNetworkPlot')
        CmdRunner.cmd(['cp -r {} {}'.format(self.output + '/step2_module_plots/vlnplot', bmk_2)], use_qsub=True)

        cmd_list = []
        for file in (Path(self.output) / 'step2_module_plots/ModuleNetworkPlot').glob('*'):
            cmd = 'convert  -resize 1024x1024 -density 600 {} {}/ModuleNetworkPlot/{}.png'. \
                format(file, bmk_2, file.name)
            cmd_list.append(cmd)
        CmdRunner.cmd(cmd_list)

        self.cfg_dict = read_cfg(self.config)
        wr = WebReport(config=self.config, template=self.template_file, report_name=self.cfg_dict['report_name'])
        wr.run_report()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--seob_obj', type=str, help='Seurat obj')
    parser.add_argument('-o', '--output', type=str, default='hdWGCNA_results', help='output dir')
    parser.add_argument('-c', '--config', type=str, help='config')

    input_args = parser.parse_args()

    hd = hdWGCNA(seob_obj=input_args.seob_obj,
                 output=input_args.output,
                 config=input_args.config)

    make_summary(Path(__file__), status='doing')
    hd.run_hdWGCNA()
    hd.run_enrich_analysis()
    hd.run_report()
    make_summary(Path(__file__), status='done')
