import sys

sys.path.append('/share/nas1/zhangjm/workspace/MyUtils')
from myrunner import MyPath, MyRunner
from pathlib import Path
import argparse
import re


class WebReport:
    """
    基于template以及config文件生成相关报告
    """
    MK_XML = '/share/nas2/genome/bin/perl {}/report_xml.pl'.format(Path(__file__).parent)
    XML_HTML = '/share/nas2/genome/biosoft/Python/2.7.8/bin/python {}/xml2HtmlConverter.py'.format(
        Path(__file__).parent)

    def __init__(self, config, template, report_name='BMK_Report'):
        """
        - config：替换的表
        - template：报告模版
        - report_name：最终压缩的报告名称
        """
        self.config = config
        self.template = template
        self.report_name = report_name

    def config_template(self):
        with open(self.config, 'r', encoding='utf-8') as f:
            config_dict = {}
            for line in f:
                k, v = line.strip().split()
                config_dict[k] = str(v)

        # 修改后的template，用于生成最后的报告
        with open('_tmp_template.txt', 'w', encoding='utf-8') as f2:
            with open(self.template, '+r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    for i in re.findall('\${[a-zA-Z|_]+}', line):
                        line = line.replace(i, config_dict[re.sub('$|{|}', '', i).replace('$', '')])

                    f2.write('{}\n'.format(line))

    @MyRunner.cmd_wrapper(threads_num=1)
    def make_report(self):
        # 一定得Web_Report
        MyPath.mkdir('Web_Report/Template')
        cmd = '{} -i Web_Report -c {} -t _tmp_template.txt -x Web_Report/configtest_xmlconvert.xml -b'. \
            format(self.MK_XML, self.config)
        cmd += ' && {} -i Web_Report/configtest_xmlconvert.xml -o Web_Report'.format(self.XML_HTML)
        cmd += ' && rm Web_Report/configtest_xmlconvert.xml _tmp_template.txt'
        cmd += ' && tar -zcvf {}.tar.gz Web_Report'.format(self.report_name)

        return [cmd]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make a web report')
    parser.add_argument('-c', '--config', type=str)
    parser.add_argument('-t', '--template', type=str)
    input_args = parser.parse_args()

    wr = WebReport(config=input_args.config, template=input_args.template)
    wr.config_template()
    wr.make_report()
