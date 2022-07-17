import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))

from bioinfor_tools.utils import MyPath, make_summary
from make_web_report.web_report import WebReport
from bioinfor_tools.cmd_runner import CmdRunner
import argparse
import logging
import shutil

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class CellTrace:
    PYTHON = '/share/nas1/zhangjm/software/miniconda3/envs/cpdb/bin/python3'
    RSCRIPT = '/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/bin/Rscript'
    run_cell_trace = Path(__file__).parent / 'run_cell_trace.R'
    template_file = Path(__file__).parent / 'template.txt'

    def __init__(self, seob,
                 cell_type='cellType',
                 output='cell_trace_results',
                 report_name='cell_trace',
                 report_title='基于monocle的轨迹分析'):

        self.seob = seob
        self.cellType = cell_type
        self.output = output
        self.report_name = report_name
        self.report_title = report_title
        self.config = None

    @CmdRunner.cmd_wrapper(use_qsub=True)
    def run_cell_trace_monocle2(self):
        if Path(self.output, 'cell_trajectory_heatmap.png').exists():
            cmd = 'echo cell trace done, skipping ...'
        else:
            cmd = '{} {} --seurat_obj {} --cell_type {} --output {}'. \
                format(self.RSCRIPT, self.run_cell_trace, self.seob, self.cellType, self.output)

        self.config = 'config.txt'

        return [cmd]

    def run_report(self):
        """生成报告，一定得Web_Report，问就是不懂"""

        with open('config_.txt', 'w', encoding='utf-8') as f:
            f.write('Project_name {}'.format(self.report_title) + '\n')
            with open(self.config, 'r') as f1:
                for i in f1:
                    f.write(i)

        self.config = 'config_.txt'
        des = 'Web_Report/BMK_1_data'
        MyPath.mkdir(des)

        for i in Path(self.output).glob('**/*'):
            if str(i).endswith('.png') or str(i).endswith('.csv') or str(i).endswith('.xls'):
                try:
                    shutil.copy(str(i), des)
                except:
                    pass

        wr = WebReport(config=self.config, template=self.template_file, report_name=self.report_name)
        wr.config_template()
        wr.make_report()
        Path('config_.txt').unlink()


if __name__ == '__main__':
    desc = """
    Version: Version v1.0
    Contact: zhangjm <zhangjm@biomarker.com.cn>
    Program Date: 2022.04.15
    Program UpDate: -
    Description: cell trace
        - based on origin pipeline
    """

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--seurat_obj', type=str, help='Seuratobj')
    parser.add_argument('-o', '--output', type=str, default='cell_trace_results', help='output dir')
    parser.add_argument('-c', '--cell_type_col', type=str, default='cellType', help='colname of celltype')
    parser.add_argument('-r', '--report_name', type=str, default='cell_trace', help='report_name')
    parser.add_argument('-t', '--report_title', type=str, default='基于monocle的轨迹分析', help='report_title')
    input_args = parser.parse_args()

    ct = CellTrace(seob=input_args.seurat_obj,
                   cell_type=input_args.cell_type_col,
                   output=input_args.output,
                   report_name=input_args.report_name,
                   report_title=input_args.report_title)

    make_summary(Path(__file__), status='doing')
    ct.run_cell_trace_monocle2()
    ct.run_report()
    make_summary(Path(__file__), status='done')