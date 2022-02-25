#!/usr/bin/env python
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent))
from myrunner import MyRunner, MyPath
from collections import defaultdict
import argparse
import logging

FORMAT = '%(asctime)s %(threadName)s -> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


@MyRunner.cmd_wrapper(threads_num=20)
def move_to_sample_dir(merge_data_path, anno_df):
    """基于样本名，将合并后的fastq移动至相应的样本名文件夹，并修改fastq名称"""
    sample_dict = {}
    merge_data = list(Path(merge_data_path).glob('*.gz'))

    with open(anno_df, 'r') as f:
        for line in f:
            try:
                sample_name, pattern = line.strip().split('\t')
                sample_dict[sample_name.strip()] = \
                    [i for i in merge_data if i.match('*{}*'.format(pattern.strip()))]
            except Exception as e:
                logging.info(e)
                logging.info('wrong anno_df file!!!Maybe in the last line.')

    cmd_list = []
    for sample_name, file_list in sample_dict.items():
        des = merge_data_path + '/' + sample_name
        MyPath.mkdir(des)
        prob = []
        for file in file_list:

            if str(file.name).endswith('_good_1.fq.gz'):
                des_name = str(file.name).replace('_good_1.fq.gz', '_S01_L001_R1_001.fastq.gz')
            elif str(file.name).endswith('_good_2.fq.gz'):
                des_name = str(file.name).replace('_good_2.fq.gz', '_S01_L001_R2_001.fastq.gz')
            elif str(file.name).endswith('.R1.fastq.gz'):
                des_name = str(file).replace('.R1.fastq.gz', '_R1_001.fastq.gz')
            elif str(file.name).endswith('.R2.fastq.gz'):
                des_name = str(file.name).replace('.R2.fastq.gz', '_R2_001.fastq.gz')
            else:
                des_name = file.name

            cmd = 'mv {} {}'.format(file, des)

            if des_name != file.name:
                rename_cmd = ' && mv {}/{} {}/{}'.format(des, file.name, des, Path(des_name).name)
                cmd += rename_cmd

            cmd_list.append(cmd)
            # 防止同一个样本的多个文件名字一样从而导致移动被覆盖
            if des_name in prob:
                logging.info('BUG: the same file name = {}'.format(file.name))
            prob.append(des_name)

    return cmd_list


@MyRunner.cmd_wrapper(threads_num=20)
def merge_fq(raw_data_path: str, merge_data_path: str):
    """将下机数据中相同的fastq文件合并在一起，输出到merge_data_path"""
    if not Path(merge_data_path).exists():
        Path(merge_data_path).mkdir()

    my_dict = defaultdict(list)
    cmd_list = []

    for i in Path(raw_data_path).glob('**/*.gz'):
        my_dict[i.name].append(str(i))

    for k, v in my_dict.items():
        if len(v) > 1:
            cmd_list.append('cat {} > {}/{}'.format(' '.join(v), merge_data_path, k))
        else:
            cmd_list.append('cp {} {}/{}'.format(v[0], merge_data_path, k))

    return cmd_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merging fastq')
    parser.add_argument('-i', '--input', type=str, help='input dir')
    parser.add_argument('-o', '--output', type=str, help='output dir', default='Merged_Data')
    parser.add_argument('-a', '--anno_df', type=str, help='sample anno ps:sample1 pattern')

    args = parser.parse_args()

    merge_fq(args.input, args.output)
    move_to_sample_dir(args.output, args.anno_df)
