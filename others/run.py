# from concurrent.futures import ThreadPoolExecutor, wait, ALL_COMPLETED, as_completed
# import subprocess
# import functools
# import time
# import logging
#
# FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
# logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
#
#
# def run_shell(cmd, test=False):
#     logging.info('subprocess is running cmd={}'.format(cmd))
#     if not test:
#         p = subprocess.Popen(cmd, shell=True)
#         ret = p.wait()
#         if ret != 0:
#             raise Exception('cmd=`{}` failed:\n{}'.format(cmd, ret))
#
#
# cmd = 'Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/TCR.R -i t_results/input/ -o t_results/TRA -t TRA -s Homo_sapiens && Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/TCR.R -i t_results/input/ -o t_results/TRB -t TRB -s Homo_sapiens && Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/Stat.R -i t_results/input/ -o t_results/stat -t TCR && export PATH=/share/nas2/genome/biosoft/gcc/5.5.0/bin/:$PATH && export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/gcc/5.5.0/lib:/share/nas2/genome/biosoft/gcc/5.5.0/lib64:/share/nas2/genome/biosoft/gsl/2.2.1/lib/:$LD_LIBRARY_PATH && Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/scRepertoire.R -i t_results/input -o t_results/combined -r  -t T'
#
# cmd2 = 'Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/BCR.R -i b_results/input/ -o b_results/IGH -t IGH -s Homo_sapiens && Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/BCR.R -i b_results/input/ -o b_results/IGK -t IGK -s Homo_sapiens && Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/BCR.R -i b_results/input/ -o b_results/IGL -t IGL -s Homo_sapiens && Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/Stat.R -i b_results/input/ -o b_results/stat -t BCR && export PATH=/share/nas2/genome/biosoft/gcc/5.5.0/bin/:$PATH && export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/gcc/5.5.0/lib:/share/nas2/genome/biosoft/gcc/5.5.0/lib64:/share/nas2/genome/biosoft/gsl/2.2.1/lib/:$LD_LIBRARY_PATH && Rscript /share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/yun_vdj/VDJ/scRepertoire.R -i b_results/input -o b_results/combined -r  -t B'
#
# # run_shell(cmd)
#
# # with ThreadPoolExecutor(max_workers=2) as executor_cmd:
# #     run_shell = functools.partial(run_shell)
# #     task_list = []
# #     for cmd in [cmd]:
# #         task = executor_cmd.submit(run_shell, cmd)
# #         task_list.append(task)
# #
# #     for future in as_completed(task_list):
# #         data = future.result()
# #
# # logging.info('-------------->')
#
# # with ThreadPoolExecutor(max_workers=2) as executor_cmd:
# #     run_shell = functools.partial(run_shell)
# #     task_list = []
# #     for cmd in [cmd2]:
# #         task = executor_cmd.submit(run_shell, cmd)
# #         task_list.append(task)
# #
# #     for future in as_completed(task_list):
# #         data = future.result()
#
# print('qc/G1Po-b_S1_L001_R2_001_fastqc,qc/G1Po-sc_S1_L001_R2_001_fastqc,qc/G1Po-t_S1_L001_R2_001_fastqc,'
#       'qc/G1Pr-b_S1_L001_R2_001_fastqc,qc/G1Pr-sc_S1_L001_R2_001_fastqc,qc/G1Pr-t_S1_L001_R2_001_fastqc,'
#       'qc/G3Po-b_S1_L001_R2_001_fastqc,qc/G3Po-sc_S1_L001_R2_001_fastqc,qc/G3Po-t_S1_L001_R2_001_fastqc,'
#       'qc/G3Pr-b_S1_L001_R2_001_fastqc,qc/G3Pr-sc_S1_L001_R2_001_fastqc,qc/G3Pr-t_S1_L001_R2_001_fastqc')
# print('cellranger/G1Po-b.origin_results,cellranger/G1Po-sc.origin_results,cellranger/G1Po-t.origin_results,'
#       'cellranger/G1Pr-b.origin_results,cellranger/G1Pr-sc.origin_results,cellranger/G1Pr-t.origin_results,'
#       'cellranger/G3Po-b.origin_results,cellranger/G3Po-sc.origin_results,cellranger/G3Po-t.origin_results,'
#       'cellranger/G3Pr-b.origin_results,cellranger/G3Pr-sc.origin_results,cellranger/G3Pr-t.origin_results')
# from pathlib import Path
# cell_cycle = 'None'
# cell_cycle_all_files = cell_cycle.split(',')
# cell_cycle = ','.join([i for i in cell_cycle_all_files if Path(i).name.startswith('analysed_integrated')])
# print(cell_cycle_all_files)
# print(cell_cycle)
a = 'a'
if a in ['a','b'] or not a == 'a':
    print('done')