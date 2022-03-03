#!/usr/bin/env python
from myrunner import MyPath, MyRunner
from sc_utils import ScBasic
from pathlib import Path
import argparse
import shutil
import re


class WebReport(ScBasic):
    produce_template_runner = Path(__file__).parent / 'WebReport/produce_template.pl'
    report_xml_runner = Path(__file__).parent / 'WebReport/report_xml.pl'
    xml2HtmlConverter_runner = Path(__file__).parent / 'WebReport/xml2HtmlConverter.py'
    package_runner = Path(__file__).parent / 'WebReport/package.pl'
    qc_report_runner = Path(__file__).parent / 'WebReport/qc_report.pl'

    dir_1 = ['BMK_1_rawData', 'BMK_2_cellranger_analysis', 'BMK_3_seurat_analysis', 'BMK_4_vdj_analysis']

    # 单细胞样本个数
    sc_sample_num = 0

    def __init__(self, config):
        super().__init__(config)
        MyPath.mkdir(*self.dir_1)

    @staticmethod
    @MyRunner.count_running_time
    def bmk_1_rawdata(gc, qc: str):
        """
        - gc: data_access_stat -> AllSample_GC_Q.stat
        - gc: data_access_fastqc/*fastqc
          - per_base_quality.png, per_sequence_gc_content.png, per_base_sequence_content.png
        """
        anchor = {'per_base_quality.png', 'per_sequence_gc_content.png', 'per_base_sequence_content.png'}
        shutil.copy(gc, 'BMK_1_rawData/{}'.format('AllSample_GC_Q.xls'))
        for orig in qc.split(','):
            sample_name = Path(orig).name.split('_')[0]
            for file in Path(orig).glob('Images/**/*'):
                if file.name in anchor:
                    shutil.copy(file, 'BMK_1_rawData/{}.{}'.format(sample_name, file.name))

    def bmk_2_cellranger_analysis(self, cellranger_sample_res: str, cellranger_stat):
        """
        所有样本的cellranger结果，包括单细胞，t细胞，b细胞
        - cellranger_sample_res -> cellranger结果，需要修改cellragner.pl将输出结果改为sample-[t|b|sc].origin_results
        """
        dir_sc = 'BMK_1_summary', 'BMK_2_PCA', 'BMK_3_TSNE', 'BMK_4_UMAP', 'BMK_5_cluster', 'BMK_6_diff'
        dir_tb = 'BMK_1_summary', 'BMK_2_clonotypes'
        des_dir_sc = ['BMK_2_cellranger_analysis/BMK_1_SC/' + i for i in dir_sc]
        des_dir_t = ['BMK_2_cellranger_analysis/BMK_2_T/' + i for i in dir_tb]
        des_dir_b = ['BMK_2_cellranger_analysis/BMK_2_B/' + i for i in dir_tb]
        MyPath.mkdir(*des_dir_sc + des_dir_t + des_dir_b)
        sample_list = cellranger_sample_res.split(',')

        for file in list(Path(cellranger_stat).glob('**/*')):
            if str(file.name).startswith('sc_'):
                shutil.copy(file, './BMK_2_cellranger_analysis/BMK_1_SC/BMK_1_summary')
            if str(file.name).startswith('b_'):
                shutil.copy(file, './BMK_2_cellranger_analysis/BMK_2_B/BMK_1_summary')
            if str(file.name).startswith('t_'):
                shutil.copy(file, './BMK_2_cellranger_analysis/BMK_2_T/BMK_1_summary')

        config_dict = {'tsne': 'BMK_3_TSNE',
                       'umap': 'BMK_4_UMAP',
                       'pca': 'BMK_2_PCA',
                       'clustering': 'BMK_5_cluster',
                       'diffexp': 'BMK_6_diff',
                       't': 'BMK_2_T',
                       'b': 'BMK_2_B'}

        for sample in sample_list:
            sample_name = Path(sample).name.replace('.origin_results', '')
            if sample_name.endswith('.sc'):
                self.sc_sample_num += 1
                shutil.copy('{}/outs/web_summary.html'.format(sample),
                            'BMK_2_cellranger_analysis/BMK_1_SC/BMK_1_summary/{}.web_summary.html'.
                            format(sample_name))

                all_file = [file for file in Path(sample).glob('outs/analysis/**/*') if file.is_file()]

                for file in all_file:
                    if file.name not in ['projection.csv', 'clusters.csv', 'differential_expression.csv']:
                        continue

                    file_name = file.parts[-2].replace('2_components', '').replace('10_components', ''). \
                                    replace('_clusters', '.').replace('_', ''). \
                                    replace('graphclust', 'graphclust.') + file.name.replace('csv', 'xls')

                    shutil.copy(file, 'BMK_2_cellranger_analysis/BMK_1_SC/{}/{}.{}'.
                                format(config_dict[file.parts[-3]], sample_name, file_name))

            else:
                des_tb = config_dict['t'] if sample_name.endswith('t') else config_dict['b']

                shutil.copy('{}/outs/web_summary.html'.format(sample),
                            'BMK_2_cellranger_analysis/{}/BMK_1_summary/{}.web_summary.xls'.
                            format(des_tb, sample_name))

                tmp_a = 'BMK_2_cellranger_analysis/{}/BMK_2_clonotypes/{}.clonotypes.xls'. \
                    format(des_tb, sample_name)
                shutil.copy('{}/outs/clonotypes.csv'.format(sample), tmp_a)

                tmp_b = 'BMK_2_cellranger_analysis/{}/BMK_2_clonotypes/{}.consensus_annotations.xls'. \
                    format(des_tb, sample_name)
                shutil.copy('{}/outs/consensus_annotations.csv'.format(sample), tmp_b)

                tmp_c = 'BMK_2_cellranger_analysis/{}/BMK_2_clonotypes/{}.filtered_contig_annotations.xls'. \
                    format(des_tb, sample_name)
                shutil.copy('{}/outs/filtered_contig_annotations.csv'.format(sample), tmp_c)

                cmd_list = ' && '.join(["sed -i 's/,/\t/g' {}".format(i) for i in (tmp_a, tmp_b, tmp_c)])

                MyRunner.runner(cmd_list=[cmd_list])

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=1)
    def bmk_3_seurat_analysis(self,
                              sample_filter,
                              sample_analysis,
                              single_enrichment,
                              single_sample_ppi,
                              single_tf_analysis,
                              single_typeanno,
                              integrated_result,
                              allcluster_statistic,
                              groupdiff_statistic,
                              cell_cycle,
                              allsample_trace,
                              cfg1,
                              cfg2):

        wang_xi_perl = Path(__file__).parent / 'WebReport/analysis_report.pl'

        # 传进来其实是所有的类型的文件，没有区分整合数据还是单个样本数据，因此在后续进行拆分并分别传参给王茜的脚本
        all_anno_file = single_enrichment.split(',')
        all_ppi_file = single_sample_ppi.split(',')
        all_tf_file = single_tf_analysis.split(',')
        all_type_anno_file = single_typeanno.split(',')

        ######### 功能富集
        single_enrichment = ','.join(
            [i for i in all_anno_file if not Path(i).name.startswith('cluster_diff_integrated') and
             not str(Path(i).name).startswith('analysed_integrated') and
             not re.findall('.*_vs_.*', str(Path(i).name))])

        # cluster_enrichment, group_enrichment
        cluster_enrichment = [i for i in all_anno_file if Path(i).name.startswith('cluster_diff_integrated')]
        if not cluster_enrichment:
            cluster_enrichment = 'None'
        else:
            cluster_enrichment = cluster_enrichment[0]

        group_enrichment = [i for i in all_anno_file if re.findall('.*_vs_.*', str(Path(i).name))]
        if not group_enrichment:
            group_enrichment = 'None'
        else:
            group_enrichment = ','.join(group_enrichment)

        ######### PPI分析
        single_sample_ppi = ','.join(
            [i for i in all_ppi_file if not Path(i).name.startswith('cluster_diff_integrated') and
             not Path(i).name.startswith('analysed_integrated') and
             not re.findall('.*_vs_.*', str(Path(i).name))])

        cluster_ppi = [i for i in all_ppi_file if Path(i).name.startswith('cluster_diff_integrated')]
        if not cluster_ppi:
            cluster_ppi = 'None'
        else:
            cluster_ppi = cluster_ppi[0]

        ######## TF分析
        single_tf_analysis = ','.join(
            [i for i in all_tf_file if not Path(i).name.startswith('cluster_diff_integrated') and
             not re.findall('.*_vs_.*', str(Path(i).name)) and
             not Path(i).name.startswith('analysed_integrated')])

        cluster_tf_result = [i for i in all_tf_file if Path(i).name.startswith('cluster_diff_integrated')]
        if not cluster_tf_result:
            cluster_tf_result = 'None'
        else:
            cluster_tf_result = cluster_tf_result[0]

        ######## 细胞注释
        single_typeanno = ','.join(
            [i for i in all_type_anno_file if not Path(i).name.startswith('analysed_integrated')])

        cell_typeanno = [i for i in all_type_anno_file if Path(i).name.startswith('analysed_integrated')]
        if not cell_typeanno:
            cell_typeanno = 'None'
        else:
            cell_typeanno = cell_typeanno[0]

        ######## 细胞周期
        if cell_cycle != 'None':
            cell_cycle_all_files = cell_cycle.split(',')
            cell_cycle = [i for i in cell_cycle_all_files if Path(i).name.startswith('analysed_integrated')]
            if not cell_cycle:
                cell_cycle = 'None'
            else:
                cell_cycle = cell_cycle[0]

        ##### 这里只拷贝多样本的
        if allsample_trace != 'None':
            cell_trace_all_files = allsample_trace.split(',')
            allsample_trace = [i for i in cell_trace_all_files if Path(i).name.startswith('analysed_integrated')]
            if not allsample_trace:
                allsample_trace = 'None'
            else:
                allsample_trace = allsample_trace[0]

        cmd = '{} {} --sample_filter {} ' \
              '--sample_analysis {} ' \
              '--single_enrichment {} ' \
              '--single_sample_ppi {} ' \
              '--single_tf_analysis {} ' \
              '--single_typeanno {} ' \
              '--integrated_result {} ' \
              '--allcluster_statistic {} ' \
              '--cluster_enrichment {} ' \
              '--cluster_ppi {} ' \
              '--cluster_tf_result {} ' \
              '--groupdiff_statistic {} ' \
              '--group_enrichment {} ' \
              '--cell_typeanno {} ' \
              '--cell_cycle {} ' \
              '--allsample_trace {} ' \
              '--cfg1 {} ' \
              '--cfg2 {} '. \
            format(self.PERL, wang_xi_perl, sample_filter,
                   sample_analysis, single_enrichment,
                   single_sample_ppi, single_tf_analysis,
                   single_typeanno, integrated_result,
                   allcluster_statistic, cluster_enrichment,
                   cluster_ppi, cluster_tf_result, groupdiff_statistic,
                   group_enrichment, cell_typeanno, cell_cycle, allsample_trace, cfg1, cfg2)

        return [cmd]

    @staticmethod
    @MyRunner.count_running_time
    def bmk_4_vdj_analysis(b_results, t_results):
        """
        - BMK_4_vdj_analysis
            - BMK_1_Tcell
                - BMK_1_TRA
                    - BMK_1_CDR3
                    - BMK_2_VDJ
                - BMK_2_TRB
                    - BMK_1_CDR3
                    - BMK_2_VDJ
                - BMK_3_AllSamples
                - TCR_*
            - BMK_2_Bcell
                - BMK_1_IGH
                    - BMK_1_CDR3
                    - BMK_2_VDJ
                - BMK_2_IGK
                    - BMK_1_CDR3
                    - BMK_2_VDJ
                - BMK_3_IGL
                    - BMK_1_CDR3
                    - BMK_2_VDJ
                - BMK_4_AllSamples
                - BCR_*
        """
        des = 'BMK_4_vdj_analysis/BMK_1_Tcell/'
        des_TRA = des + 'BMK_1_TRA'
        des_TRB = des + 'BMK_1_TRB'
        des_all_sample = des + 'BMK_3_AllSamples'
        MyPath.mkdir(des_all_sample)
        MyPath.mkdir(*['{}/{}'.format(i, _) for i in [des_TRA, des_TRB] for _ in ('BMK_1_CDR3', 'BMK_2_VDJ')])

        des2 = 'BMK_4_vdj_analysis/BMK_2_Bcell/'
        des2_BMK_1_IGH = des2 + 'BMK_1_IGH'
        des2_BMK_2_IGK = des2 + 'BMK_2_IGK'
        des2_BMK_3_IGL = des2 + 'BMK_3_IGL'
        des2_BMK_4_AllSamples = des2 + 'BMK_4_AllSamples'
        MyPath.mkdir(des2_BMK_4_AllSamples)
        MyPath.mkdir(*['{}/{}'.format(i, _) for i in [des2_BMK_1_IGH, des2_BMK_2_IGK, des2_BMK_3_IGL] for _ in
                       ('BMK_1_CDR3', 'BMK_2_VDJ')])

        for file in list(Path(t_results).glob('**/*')) + list(Path(b_results).glob('**/*')):
            file_name = str(file.name)

            # T-cell
            if file_name.startswith('TRA_CDR3_lens') or file_name.startswith('TRA_Clonotype_Relative_abundance') or \
                    file_name.startswith('TRA_diversity_measure') or file_name.startswith('TRA_motif_') or \
                    file_name.startswith('TRA_trackClonotypes'):
                shutil.copy(file, '{}/BMK_1_CDR3'.format(des_TRA))

            if file_name.startswith('TRB_CDR3_lens') or file_name.startswith('TRB_Clonotype_Relative_abundance') or \
                    file_name.startswith('TRB_diversity_measure') or file_name.startswith('TRB_motif_') or \
                    file_name.startswith('TRB_trackClonotypes'):
                shutil.copy(file, '{}/BMK_1_CDR3'.format(des_TRB))

            if file_name.startswith('TRA_CDR3_J_') or file_name.startswith('TRA_CDR3_V_') or \
                    file_name.startswith('TRA_Circos_') or file_name.startswith('TRA_TRAJ_') or \
                    file_name.startswith('TRA_V_J_pair_'):
                shutil.copy(file, '{}/BMK_2_VDJ'.format(des_TRA))

            if file_name.startswith('TRB_CDR3_J_') or file_name.startswith('TRB_CDR3_V_') or \
                    file_name.startswith('TRB_Circos_') or file_name.startswith('TRB_TRBJ_') or \
                    file_name.startswith('TRB_V_J_pair_'):
                shutil.copy(file, '{}/BMK_2_VDJ'.format(des_TRB))

            if file_name.startswith('TCR_cloneType_tsne') or file_name.startswith('TCR_cloneType_ump') or \
                    file_name.startswith('TCR_cluster_cloneType') or file_name.startswith('TCR_cluster_Overlap'):
                shutil.copy(file, des)

            if file_name.startswith('TCR.All_Sample_diversity_measure') or \
                    file_name.startswith('TCR.All_Sample_Rank-abundance') or \
                    file_name.startswith('TCR.All_Sample_Rank-abundance') or \
                    file_name.startswith('TCR.All_Sample_rarefaction') or \
                    file_name.startswith('TCR.All_Samples_stat'):
                shutil.copy(file, des_all_sample)

            # B-cell
            if file_name.startswith('IGH_CDR3_lens') or file_name.startswith('IGH_Clonotype_Relative_abundance') or \
                    file_name.startswith('IGH_diversity_measure') or file_name.startswith('IGH_motif_') or \
                    file_name.startswith('IGH_trackClonotypes'):
                shutil.copy(file, '{}/BMK_1_CDR3'.format(des2_BMK_1_IGH))

            if file_name.startswith('IGK_CDR3_J_') or file_name.startswith('IGK_CDR3_V_') or \
                    file_name.startswith('IGK_Circos_') or file_name.startswith('IGK_IGKJ_') or \
                    file_name.startswith('IGK_V_J_pair_'):
                shutil.copy(file, '{}/BMK_1_CDR3'.format(des2_BMK_2_IGK))

            if file_name.startswith('IGL_CDR3_J_') or file_name.startswith('IGL_CDR3_V_') or \
                    file_name.startswith('IGL_Circos_') or file_name.startswith('IGL_IGLJ_') or \
                    file_name.startswith('IGL_V_J_pair_'):
                shutil.copy(file, '{}/BMK_1_CDR3'.format(des2_BMK_3_IGL))

            if file_name.startswith('IGH_CDR3_J_') or file_name.startswith('IGH_CDR3_V_') or \
                    file_name.startswith('IGH_Circos_') or file_name.startswith('IGH_IGHJ_') or \
                    file_name.startswith('IGH_V_J_pair_'):
                shutil.copy(file, '{}/BMK_2_VDJ'.format(des2_BMK_1_IGH))

            if file_name.startswith('IGK_CDR3_J_') or file_name.startswith('IGK_CDR3_V_') or \
                    file_name.startswith('IGK_Circos_') or file_name.startswith('IGL_IGLJ_') or \
                    file_name.startswith('IGK_V_J_pair_'):
                shutil.copy(file, '{}/BMK_2_VDJ'.format(des2_BMK_2_IGK))

            if file_name.startswith('IGL_CDR3_J_') or file_name.startswith('IGL_CDR3_V_') or \
                    file_name.startswith('IGL_Circos_') or file_name.startswith('IGL_IGLJ_') or \
                    file_name.startswith('IGL_V_J_pair_'):
                shutil.copy(file, '{}/BMK_2_VDJ'.format(des2_BMK_3_IGL))

            if file_name.startswith('BCR_cloneType_tsne') or file_name.startswith('BCR_cloneType_ump') or \
                    file_name.startswith('BCR_cluster_cloneType') or file_name.startswith('BCR_cluster_Overlap'):
                shutil.copy(file, des2)

            if file_name.startswith('BCR.All_Sample_diversity_measure') or \
                    file_name.startswith('BCR.All_Sample_Rank-abundance') or \
                    file_name.startswith('BCR.All_Sample_Rank-abundance') or \
                    file_name.startswith('BCR.All_Sample_rarefaction') or \
                    file_name.startswith('BCR.All_Samples_stat'):
                shutil.copy(file, des2_BMK_4_AllSamples)

    @MyRunner.count_running_time
    def end(self):
        cmd = 'mkdir Web_Report && mv BMK_* Web_Report'
        MyRunner.runner([cmd])

        # 第一步构建template
        cmd = '{} {} -i Web_Report -o Web_Report -c {} -t '. \
            format(self.PERL,
                   self.produce_template_runner,
                   self._config)
        if self.sc_sample_num > 1:
            cmd += 'm'
        else:
            cmd += 's'

        # 第二步
        cmd += ' && {} {} -i Web_Report -c {} -t ' \
               'Web_Report/template -x Web_Report/configtest_xmlconvert.xml -b'. \
            format(self.PERL,
                   self.report_xml_runner,
                   self._config)

        # 第三步
        cmd += ' && {} {} -i Web_Report/configtest_xmlconvert.xml -o Web_Report'. \
            format(self.python2, self.xml2HtmlConverter_runner)

        # 第四步
        cmd += ' && {} {} -in Web_Report -od Needed_Data -cfg {}'. \
            format(self.PERL, self.package_runner, self._config)

        MyRunner.runner([cmd])

    @MyRunner.cmd_wrapper(threads_num=1)
    def qc_report(self):
        cmd = '{} {} --input Web_Report --output {} --cfg {} --samples_num {}'. \
            format(self.PERL,
                   self.qc_report_runner,
                   './',
                   self._config,
                   self.sc_sample_num)
        cmd += '&& {} {} -i Web_Report/configtest_qc.xml -o Web_Report/QC_HTML'. \
            format(self.python2, self.xml2HtmlConverter_runner)

        return [cmd]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='running sb web report')
    parser.add_argument('--gc', type=str)
    parser.add_argument('--qc', type=str)
    parser.add_argument('--cellranger', type=str)
    parser.add_argument('--cellranger_stat', type=str)
    parser.add_argument('--sample_filter', type=str)
    parser.add_argument('--sample_analysis', type=str)
    parser.add_argument('--single_enrichment', type=str)
    parser.add_argument('--single_sample_ppi', type=str)
    parser.add_argument('--single_tf_analysis', type=str)
    parser.add_argument('--single_typeanno', type=str)
    parser.add_argument('--integrated_result', type=str, default='None')
    parser.add_argument('--allcluster_statistic', type=str, default='None')
    parser.add_argument('--groupdiff_statistic', type=str, default='None')
    parser.add_argument('--cell_cycle', type=str, default='None')
    parser.add_argument('--allsample_trace', type=str, default='None')
    parser.add_argument('--cfg1', type=str)  # data.cfg
    parser.add_argument('--cfg2', type=str)  # detail.cfg
    parser.add_argument('--t_results', type=str, default='None')
    parser.add_argument('--b_results', type=str, default='None')
    input_args = parser.parse_args()

    web = WebReport(config=input_args.cfg2)

    web.bmk_1_rawdata(gc=input_args.gc, qc=input_args.qc)
    web.bmk_2_cellranger_analysis(cellranger_sample_res=input_args.cellranger,
                                  cellranger_stat=input_args.cellranger_stat)
    web.bmk_3_seurat_analysis(sample_filter=input_args.sample_filter,
                              sample_analysis=input_args.sample_analysis,
                              single_enrichment=input_args.single_enrichment,
                              single_sample_ppi=input_args.single_sample_ppi,
                              single_tf_analysis=input_args.single_tf_analysis,
                              single_typeanno=input_args.single_typeanno,
                              integrated_result=input_args.integrated_result,
                              allcluster_statistic=input_args.allcluster_statistic,
                              groupdiff_statistic=input_args.groupdiff_statistic,
                              cell_cycle=input_args.cell_cycle,
                              allsample_trace=input_args.allsample_trace,
                              cfg1=input_args.cfg1,
                              cfg2=input_args.cfg2)
    web.bmk_4_vdj_analysis(t_results=input_args.t_results,
                           b_results=input_args.b_results)
    web.end()
    web.qc_report()
