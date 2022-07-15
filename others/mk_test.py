# import re
#
# from myrunner import MyPath, MyRunner
from pathlib import Path
import scanpy as sc
sc.read_visium()
#
#
# @MyRunner.count_running_time
# @MyRunner.cmd_wrapper()
# def mk_test_data(dir_):
#     file_dir = [Path(dir_).parent / str(i).replace('real_data', 'test_data') for i in Path(dir_).glob('*')]
#     MyPath.mkdir(*file_dir)
#     cmd_list = []
#     for file in Path(dir_).glob('**/*'):
#         if file.is_file():
#             des = str(file).replace('real_data', 'test_data')
#             cmd = 'zcat {} | head -50000000 | gzip - > {}'.format(file, des)
#             cmd_list.append(cmd)
#
#     return cmd_list
#
#
# a = '/share/nas1/zhangjm/workspace/PipeLine/test_2021.11.09/real_data'
# # mk_test_data(a)
#
# print(re.findall('.*_b_.*', str(Path('G1Po-sc_vs_G3Pr-sc.ppi_result'))))

print('asdfa&asdfad&_vs_asdf&asdfa'.replace('&', '_'))
print(Path('abd').exists())