#!/usr/bin/env python

from collections import defaultdict
from pathlib import Path
from myrunner import *
import argparse
import math



FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class FastqFilter:
    samtools = '/share/nas1/zhangjm/software/miniconda3/bin/samtools'
    bamToFastq = '/share/nas2/genome/bin/bamToFastq'
    seqkit = '/share/nas1/wuhm/bin/seqkit'

    def __init__(self, origin_bam, origin_fastq_dir, output, mapping_rate):
        self.origin_bam = origin_bam
        self.origin_fastq_dict = self.origin_fastq(origin_fastq_dir)
        self.output = output
        self.mapping_rate = mapping_rate.split(',')

    @staticmethod
    def origin_fastq(origin_fastq_dir):
        file_dict = defaultdict(list)
        for fastq in list(Path(origin_fastq_dir).glob('**/*')):
            if fastq.match('*_R1_*'):
                file_dict['R1'].append(str(fastq))
            elif fastq.match('*_R2_*'):
                file_dict['R2'].append(str(fastq))
            else:
                continue

        return file_dict

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=1)
    def get_unmapped_fastq(self):
        MyPath.mkdir(self.output)
        if (Path(self.output) / 'unmapped.R2.fastq').exists():
            cmd = 'echo {} exists'.format((Path(self.output) / 'unmapped.R2.fastq'))

        else:
            cmd = '{} view -@ 15 -b -f 4 {} > {}/file_unmapped.bam'. \
                format(self.samtools, self.origin_bam, self.output)
            cmd += '&& {} -i {}/file_unmapped.bam -fq {}/unmapped.R1.fastq -fq2 {}/unmapped.R2.fastq'. \
                format(self.bamToFastq, self.output, self.output, self.output, self.output)

        return [cmd]

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=1)
    def get_filtered_id(self):
        if Path('{}/awk_results.txt'.format(self.output)).exists():
            cmd2 = 'echo {} exists'.format((Path(self.output) / 'awk_results.txt'))

        else:
            origin_rate, des_rate = [float(i) for i in self.mapping_rate]
            logging.info('rate from {} to {}'.format(origin_rate, des_rate))

            cmd = "grep '^@' {}/unmapped.R2.fastq | wc  -l ".format(self.output)
            unmapped_reads_num = int(subprocess.check_output(cmd, shell=True).strip())
            # unmapped_reads_num = 100000

            logging.info('unmapped_reads_num={}, counting reads number to filtering'.format(unmapped_reads_num))

            total = unmapped_reads_num / (1 - origin_rate)
            mapped_reads_num = total * origin_rate
            filter_num = math.floor(total - (mapped_reads_num / des_rate))

            logging.info('total reads={};filter_num={},get filtered reads id ...'.format(total, filter_num))

            cmd2 = "grep '^@' {}/unmapped.R2.fastq | sed -e 's/@//g' -e '{}' | shuf -n{} > {}/filter_list.txt". \
                format(self.output, "s/\/2//g", filter_num, self.output)

            cmd2 += " && zcat {} | grep '^@' | sed -e 's/@//g' | cut -d' ' -f 1 > {}/total_list.txt". \
                format(' '.join(self.origin_fastq_dict['R2']), self.output)

            cmd2 += " && awk '{}' {}/filter_list.txt {}/total_list.txt > {}/awk_results.txt". \
                format('NR==FNR{ a[$1]=$1 } NR>FNR{ if(a[$1] == ""){ print $1}}',
                       self.output, self.output, self.output)

        return [cmd2]

    @MyRunner.count_running_time
    @MyRunner.cmd_wrapper(threads_num=2)
    def filter_fastq(self):
        if list(Path(self.output).glob('*.gz')):
            cmd_list = ['echo fastq files exist']

        else:
            cmd_list = []
            for _ in self.origin_fastq_dict.values():
                for file in _:
                    cmd = '{} grep --pattern-file {}/awk_results.txt --threads 10 {} -o {}/{}'. \
                        format(self.seqkit, self.output, file, self.output, Path(file).name)
                    cmd_list.append(cmd)

        return cmd_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter fastq')
    parser.add_argument('-b', '--origin_bam', type=str, help='input origin_bam')
    parser.add_argument('-f', '--origin_fastq_dir', type=str, help='origin_fastq_dir')
    parser.add_argument('-o', '--output', type=str, help='output dir')
    parser.add_argument('-m', '--mapping_rate', type=str, help='mapping_rate [0.6,0.8]')
    args = parser.parse_args()

    fq = FastqFilter(origin_bam=args.origin_bam,
                     origin_fastq_dir=args.origin_fastq_dir,
                     output=args.output,
                     mapping_rate=args.mapping_rate)

    fq.get_unmapped_fastq()
    fq.get_filtered_id()
    fq.filter_fastq()
