import sys
from pathlib import Path

sys.path.append(str(Path(__file__).absolute().parent.parent))
from bioinfor_tools.cmd_runner import CmdRunner
from bioinfor_tools.utils import MyPath
from collections import defaultdict
import pandas as pd
import argparse
import logging


class EnrichAnalysis:
    _ = Path(__file__).parent
    # RSCRIPT = '/share/nas1/zhangjm/software/miniconda3/envs/RNA_velocyto/bin/Rscript'
    RSCRIPT = 'Rscript'
    enrich_analysis = _ / 'enrich_analysis.R'

    def __init__(self,
                 diff_results_table,
                 output='GO_KEGG_results',
                 detail_config='', GO_info='', KEGG_info='',
                 title='DE_genes',
                 gene_column='gene',
                 split_by='',
                 symbol_list=''):
        self.table = diff_results_table
        self.output = output
        self.title = title
        self.gene_column = gene_column
        self.GO_info = GO_info
        self.KEGG_info = KEGG_info
        self.split_by = split_by
        self.symbol_list = symbol_list
        MyPath.mkdir(self.output)

        if not self.GO_info or not self.KEGG_info:
            self.detail_dict = self.read_cfg(detail_config)
            self.get_info()

    @staticmethod
    def read_cfg(detail_config):
        with open(detail_config, 'r', encoding='utf-8') as f:
            tmp_dict = {}
            for line in f:
                if len(line.split()) == 1:
                    tmp_dict[line.strip()] = ''
                else:
                    k, v = line.split()
                    tmp_dict[k] = v.strip()

        return tmp_dict

    def get_info(self):
        # {gene_id:[symbol1, symbol2, ...]}
        with open(self.symbol_list, 'r') as f:
            gene_id2symbol_dict = defaultdict(list)
            for _, i in enumerate(f):
                if _ == 0:
                    continue
                gene_id, symbol = i.split()
                gene_id2symbol_dict[gene_id.strip()].append(symbol.strip().replace('_', '-'))

        # 创建KEGG.info
        if not self.KEGG_info:
            kegg_info = '{}/KEGG.info'.format(self.output)
            if not Path(kegg_info).exists():
                kegg = Path(self.detail_dict['Known_anno']) / 'Known.longest_transcript.fa.Kegg.pathway'

                if not kegg.exists():
                    raise Exception('{} not exists'.format(kegg))

                with open(kegg, 'r') as f:
                    with open(kegg_info, 'w') as f1:
                        for _, i in enumerate(f):
                            if _ == 0:
                                continue
                            pathway, pathway_id, Gene_id = [i.split('\t')[_] for _ in [0, 1, 3]]
                            for gene_id in Gene_id.split(';'):
                                for symbol in gene_id2symbol_dict[gene_id]:
                                    f1.write('{}\t{}\t{}\t{}\n'.format(pathway_id, pathway, symbol, gene_id))
                logging.info('done! KEGG.info')

            self.KEGG_info = kegg_info

        # 创建GO.info
        if not self.GO_info:
            go_info = '{}/GO.info'.format(self.output)
            if not Path(go_info).exists():
                go = Path(self.detail_dict['Known_anno']) / 'Known.longest_transcript.fa.GO.anno.txt'

                if not go.exists():
                    raise Exception('{} not exists'.format(go))

                with open(go, 'r') as f:
                    with open(go_info, 'w') as f1:
                        for _, i in enumerate(f):
                            if _ == 0:
                                continue
                            gene_id, go_anno = i.strip().split('\t')[0], i.strip().split('\t')[2::]
                            for _ in go_anno:
                                go_class = _.split(':')[0]
                                go_iterm = _.split(':')[1].replace('(GO', '').strip()
                                go_id = 'GO:' + _.split(':')[2].replace(');', '').strip()
                                for symbol in gene_id2symbol_dict[gene_id]:
                                    f1.write('{}\t{}\t{}\t{}\t{}\n'.format(go_id, go_iterm, symbol, go_class, gene_id))
                logging.info('done! GO.info')

            self.GO_info = go_info

    def run_enrichment(self):
        if not self.split_by:
            logging.info('doing enrich analysis')
            cmd = '{} {} --table {} --gene_colname {} --output {} --title {}'. \
                format(self.RSCRIPT,
                       self.enrich_analysis,
                       self.table,
                       self.gene_column,
                       self.output,
                       self.title)
            if self.GO_info:
                cmd += ' --GO.info {}'.format(self.GO_info)
            if self.KEGG_info:
                cmd += ' --KEGG.info {}'.format(self.KEGG_info)

            CmdRunner.cmd([cmd], use_qsub=True)

        else:
            data = pd.read_csv(self.table, sep='\t')
            cmd_list = []
            logging.info('splitting ...')
            for cluster, group in data.groupby(self.split_by):
                cluster = str(cluster)
                for i in ' *&.,':
                    cluster = cluster.replace(i, '_')
                res = '{}/{}.genelist.xls'.format(self.output, cluster)
                group.to_csv(res, sep='\t', index=False)

                cmd = '{} {} --table {} --gene_colname {} --output {} --title {}'. \
                    format(self.RSCRIPT,
                           self.enrich_analysis,
                           res,
                           self.gene_column,
                           (Path(self.output) / '{}.results'.format(cluster)),
                           self.title)
                if self.GO_info:
                    cmd += ' --GO.info {}'.format(self.GO_info)
                if self.KEGG_info:
                    cmd += ' --KEGG.info {}'.format(self.KEGG_info)
                cmd_list.append(cmd)

            CmdRunner.cmd(cmd_list, use_qsub=True, n_jobs=10)

            # 满足sb逻辑
            # for file in Path(self.output).glob('**/*'):
            #     if file.is_file():
            #         shutil.copy(file, str(file).replace('/', '_'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-t', '--table', type=str, help='input table')
    parser.add_argument('-o', '--output', type=str, help='colnames of cell type', default='GO_KEGG_results')
    parser.add_argument('-c', '--detail_config', type=str, help='detail config')
    parser.add_argument('-g', '--GO_info', type=str, help='GO info')
    parser.add_argument('-k', '--KEGG_info', type=str, help='KEGG info')
    parser.add_argument('-i', '--title', type=str, help='title', default='DE_genes')
    parser.add_argument('-u', '--column', type=str, help='column of gene', default='gene')
    parser.add_argument('-s', '--split_by', type=str, help='column of gene', default='cluster')
    parser.add_argument('-a', '--symbol_list', type=str, help='symbol.list', default='')
    input_args = parser.parse_args()

    ea = EnrichAnalysis(diff_results_table=input_args.table,
                        output=input_args.output,
                        detail_config=input_args.detail_config,
                        GO_info=input_args.GO_info,
                        KEGG_info=input_args.KEGG_info,
                        title=input_args.title,
                        gene_column=input_args.column,
                        split_by=input_args.split_by,
                        symbol_list=input_args.symbol_list)
    ea.run_enrichment()
