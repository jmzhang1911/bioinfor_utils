#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import sys, os, argparse, glob, os.path, time, glob

reload(sys)
sys.setdefaultencoding('utf8')
import numpy as np
import random
import math
import re
import xml.dom.minidom
from xml.dom.minidom import parseString, getDOMImplementation
from collections import OrderedDict
import urllib
import shutil
import codecs
import subprocess

Bin = os.path.split(os.path.realpath(__file__))[0]
t1 = 0
t2 = 0
t3 = 0


class parseXML(object):
    def __init__(self, file):
        self.f = file

    def __iter__(self):
        return self

    def parse(self):
        file = open(self.f, 'r')

    def next(self):
        r = [self.s, self.e]
        self.s, self.e = self.bin + self.s, self.e + self.bin
        return r


def cp_template(o, s=None):
    if os.path.exists(o + '/src'): shutil.rmtree(o + '/src')
    if s:
        shutil.copytree(s, o + '/src')
    else:
        shutil.copytree(Bin + '/src', o + '/src')


def convert_pic():
    shutil.copyfile("oldfile", "newfile")


def addTable(tmp, isfirst):
    result = ""
    if isfirst:
        for i in tmp:
            result += "<th>%s</th>" % i
    else:
        for i in tmp:
            result += "<td>%s</td>" % i
    return result


def timenow():
    """return current time as a string
    """
    return time.strftime('[%d/%m/%Y %H:%M:%S]', time.localtime(time.time()))


def Indent(dom, node, indent=0):
    children = node.childNodes[:]
    if indent:
        text = dom.createTextNode('\n' + '\t' * indent)
        node.parentNode.insertBefore(text, node)
    if children:
        if children[-1].nodeType == node.ELEMENT_NODE:
            text = dom.createTextNode('\n' + '\t' * indent)
            node.appendChild(text)
        for n in children:
            if n.nodeType == node.ELEMENT_NODE:
                Indent(dom, n, indent + 1)


def addElement(d, r, l, **args):
    global t1
    global t2
    global t3
    if (re.match("^h1$", l)):
        t1 += 1
        t2 = 0
        t3 = 0
    if (re.match("^h2$", l)):
        t2 += 1
        t3 = 0
    if (re.match("^h3$", l)):
        t3 += 1
    if t2 == 0:
        myleble = "%s" % (t1)
    elif t3 == 0:
        myleble = "%s.%s" % (t1, t2)
    else:
        myleble = "%s.%s.%s" % (t1, t2, t3)

    e = d.createElement(l)
    for k, v in args.iteritems():
        if (re.match("^h[1-3]$", l) and re.match("^name$", k)):
            e.setAttribute(k, "%s %s" % (myleble, v))
        else:
            e.setAttribute(k, v)
    r.appendChild(e)
    return e


def addElementListImage(d, r, l, p, analysis_path, desc=None):
    if desc == None:
        desc = r.getAttribute('desc')
    for i in sorted(p):
        dir, suffix = os.path.splitext(i)
        name = os.path.basename(dir)
        addElement(d, r, l, name=name, type="type1", path=i.replace(analysis_path + '/', ""), desc=desc)


def get_reflen(refinfo):
    print
    refinfo
    r = re.split('\s+', refinfo[-1])[1]
    print
    r
    m = re.match('(\d*\.?\d*?)([kKmMGg]?)', r)
    if m:
        reflen = m.group(1)
        k = 1
        if not m.group(2) == '':

            if re.match('[kK]', m.group(2)): k = 1000
            if re.match('[mM]', m.group(2)): k = 1000000
            if re.match('[gG]', m.group(2)): k = 1000000000
        return float(reflen) * k
    else:
        sys.stderr.write("Can't find reference length please check your config file\n")
        exit(1)


parser = argparse.ArgumentParser(description='This script was used to generate a xml report from reseq analysis.')
parser.add_argument('-cfg', '--cfg', dest='cfg', required=True, help='input config file')
parser.add_argument('-l', '--tableLineNum', dest='tableLineNum', required=False, type=int, default=100,
                    help='Input max table line numbers to show, default:100')
parser.add_argument('-n', '--name', dest='name', required=False, default='configtest',
                    help='specify the output file prefix,default is "configtest"')
args = parser.parse_args()
###################################################################
lineNum = args.tableLineNum
###########################################################################
f = codecs.open(args.cfg, 'rU', 'utf-8')
for i in f:
    i = i.strip()
    if i.startswith('species'): species = re.split(u'\s+', i)[1]
    if i.startswith('title'): title = re.split(u'\s+', i)[1]
    if i.startswith('sampleNum'): sampleNum = re.split(u'\s+', i)[1]
    if i.startswith('total'): total = re.split(u"\s+", i)[1]
f.close()
#####################构造图片
os.system("/share/nas2/genome/biosoft/Python/2.7.8/bin/python " + Bin + "/data/draw_pic.py")

#####################################get project guideLine#############################################
cleanData = int(total) * int(sampleNum)
averagemapRate = float(np.random.uniform(90, 99, 1))
averageDepth = int(np.random.randint(10, 15, 1))
coverageRate = float(np.random.uniform(80, 95, 1))
#####################################################################################
impl = xml.dom.minidom.getDOMImplementation()
dom = impl.createDocument(None, 'report', None)
root = dom.documentElement
addElement(dom, root, 'report_version', value="v1.3")
addElement(dom, root, 'report_name', value='10X genomics单细胞转录组RNA速度分析结题报告')
addElement(dom, root, 'report_code', value=species)
addElement(dom, root, 'report_user', value="XXXX")
addElement(dom, root, 'report_user_addr', value="XXXX")
addElement(dom, root, 'report_time',
           value="XXXX/XX/XX;XXXX/XX/XX;XXXX/XX/XX;%s" % time.strftime('%Y/%m/%d', time.localtime(time.time())))

# 摘要
cell_number = sum([random.randint(50000, 51000) for i in range(int(sampleNum))])
data_size = sum([random.randint(48, 54) for i in range(int(sampleNum))])
addElement(dom, root, 'report_abstract', value="<p class=\" p-abstract\" >"
                                               "使用Velocyto软件及scVelo对%s个单细胞样本进行RNA速度分析，基于mRNA的spliced / unspliced信息，推断RNA变化速率和方向，并预测细胞的未来状态和最终命。"
                                               % ('6'))

# 1.背景介绍
des0 = """
单细胞测序分许可以解决细胞异质性的问题，特别是拟时序分析通过构建细胞间变化轨迹预测细胞随时间的变化。
以细胞的表达量数据为研究对象，采用monocle3软件在虚拟时间轴上对细胞的变化模式进行分析，模拟重建细胞的动态变化过程，
获得细胞间的状态转换关系，以及不同状态细胞间差异基因的表达情况。但是拟时序分析（伪时间分析）无法精准判断细胞分化的方向，
只针对基因表达水平相似度无法精准真实还原细胞状态和分化路径。"""
des0_1 = """
捕获细胞分化趋势的能力可有助于科学家更好地分析复杂组织和器官中的细胞功能和功能障碍，针对细胞mRNA降解速度模拟推断细胞分化的趋势，从而揭示细胞分化轨迹。
由2018年由瑞典卡罗林斯卡医学院的Sten Linnarsson和 Peter V. Kharchenko[1]实验室提出一种新概念：RNA velocity （RNA速度），
通过对unspliced(nascent) 和spliced(mature) mRNA丰度的评估在时间维度上揭示转录本动态变化的一个指标，利用动态模型评估细胞瞬间状态的RNA速度，一种动态对细胞状态的观测，可以推断细胞轨迹。
采用Velocyto(http://velocyto.org/)软件可以估算RNA速度，推断细胞分化轨迹。该分析研究意义主要在于：<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 1) 根据检测细胞的类型推断细胞分化的轨迹 <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 2) 预测单个细胞的变化方向，得到细胞间的转变过程 <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 3) 对以往数据进行挖掘，采用初生（未剪接的）和成熟（剪接的）mRNA的相对丰度用来估计基因剪接和降解的速率 <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 4) 未稳定的细胞向处于稳泰平衡的细胞分化，揭示动力学稳态研究过程 <br />
"""

addElement(dom, root, 'h1', name="背景介绍", type="type1", desc="背景介绍")
addElement(dom, root, 'p', type="type1", desc=des0)
addElement(dom, root, 'p', type="type1", desc=des0_1)

# 2. 分析方法
des1_1 = """
RNA丰度是单个细胞状态的有力指标，而在整个生命周期中，细胞中含有处于不同阶段的mRNA混合物---新生的mRNA前体、成熟的mRNA和片段化的遭受切割的mRNA，
新生的mRNA和成熟mRNA的降解一直都是趋于平衡，任何方向的转变---过多的或过少的新生mRNA被移除---都预示着细胞行为的变化。
通过Velocyto软件计算初生(未剪接的)和成熟(剪接的)mRNA的相对丰度可以用来估计基因剪接和降解的速率，在这个动态过程中，
转录速率的增加导致unspliced mRNA的快速增加, 其次是增加subsequent拼接mRNA，直到达到一个新的稳定状态。
RNA速度（RNA velocity）---RNA随时间变化的速率，这种RNA速度可在以小时计的尺度上作为细胞命运的预测因子。
"""
des1_2 = """
分析原理：Velocyto软件计算每个细胞未剪接RNA（u）和剪接RNA（s）的比例，计算mRNA相对丰度变化，
γ代表降解速率，根据时间上的ds/dt的变化值，当u>γs时，说明细胞存在大量的未剪接的mRNA，
细胞处于一个初生状态，mRNA需要转录剪接分化；当u=γs，说明细胞中未剪接RNA（u）和剪接RNA（s）处于一个平稳状态，
即细胞处于一个成熟稳定期；当u<γs，说明细胞存在大量剪接的mRNA，细胞处于一个晚期的状态，mRNA需要降解。
根据RNA速度的变化判断细胞分化轨迹的推断。具体分析步骤如下：<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 1) 剪接计数和未剪接计数是通过包含内含子序列的单独计数读数来估计的 <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 2) 转录模式的动态过程模型：转录速率(α),剪接速率(β)和降解速率(γ)，未剪接RNA分析丰度（unspliced ，u)和成熟mRNA分子丰度（spliced，s）<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 3) 根据模型函数，显示未剪接和剪接的mRNA随时间响应α的动态变化 <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 4) 根据未剪接RNA分析丰度和成熟mRNA分子丰度, 预测未剪接的mRNA下一次剪接的时间点 <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 5) 关键基因在24h内剪接和未剪接mrna的丰度变化 <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 6) 根据RNA速度模型预测的未来时间t表达状态的变化
"""
addElement(dom, root, 'h1', name="分析方法", type="type1", desc="分析方法")
addElement(dom, root, 'p', type="type1", desc=des1_1)
addElement(dom, root, "pic", name="图1 RNA速度分析原理示意图", type="type1", desc=" ", path="src/images/*.png")
#addElementListImage(dom, root, "pic", desc="图2 RNA速度分析原理示意图", analysis_path="src/images/list")
addElement(dom, root, 'p', type="type1", desc=des1_2)

## 2.1 单细胞悬液的制备和检测
des2 = """
根据组织或细胞样品的来源和特性，选取适合的方法将其在较短时间内制备成单细胞悬液，
通过细胞计数和细胞活性检测，严格把控单细胞悬液的细胞浓度和细胞活性，
以保障使用合格的样品进行10X单细胞测序。如果不需要百迈客制备单细胞悬液则跳过该步骤。
"""
addElement(dom, root, 'h2', name="单细胞悬液的制备和检测", type="type1", desc="单细胞悬液的制备和检测")
addElement(dom, root, 'p', type="type1", desc=des2)

## 2.2 单细胞文库的构建和质检
des3 = """
单细胞悬液样本对细胞活率的要求为 90%，细胞浓度在 700-1200 个/μl左右。
通过微流控技术，在芯片上，将单个细胞及反应所需试剂，与带有细胞标签序列（cell Barcode）的胶珠（bead）一起包裹在GEMs液滴中，
收集包有细胞的GEMs液滴，在GEMs液滴中，细胞裂解释放出 mRNA，mRNA 与 bead 上的 cell Barcode引物结合，
完成反转录反应，随后打破GEMs，回收并通过PCR扩增富集cDNA，进行 cDNA 的文库构建。
"""
addElement(dom, root, 'h2', name="单细胞文库的构建和质检", type="type1", desc="单细胞文库的构建和质检")
addElement(dom, root, 'p', type="type1", desc=des3)

## 2.3 文库质控
des4 = """
文库质检主要包括三个方面：

（1）Qubit 4.0 检测cDNA产物和文库浓度；

（2）Agilent 2100 检测cDNA文库的插入片段大小；

（3）文库需严格满足以下条件方可进行测序：插入片段大小合格，峰型单一，无杂峰，无接头和无引物二聚体。
"""
addElement(dom, root, 'h2', name="文库质控", type="type1", desc="文库质控")
addElement(dom, root, 'p', type="type1", desc=des4)

## 2.4 上机测序
addElement(dom, root, 'h2', name="上机测序", type="type1", desc="上机测序")
addElement(dom, root, 'p', type="type1", desc='库检合格后，使用Illumina平台进行测序。')

# 3. 信息分析流程
des5 = """
利用10X Genomics官方软件CellRanger<a href=\"#ref1\">[1]</a>进行数据比对、基因定量以及细胞鉴定。
同时cellranger还基于鉴定的细胞进行了细胞亚群分析以及差异表达基因分析。
根据cellranger鉴定出的细胞类型以及定量结果，后续采用Seurat<a href=\"#ref1\">[2]</a>软件对数据做了二次过滤，
对过滤后的数据进行二次分析，包括：细胞亚群分析、差异表达基因分析、GO/KEGG富集分析以及PPI网络分析等等。
为了保证报告的流畅度，所有文件链接以及表格都只展示前6行，完整内容请见相应的结果文件。
"""

addElement(dom, root, 'h1', name="信息分析流程", type="type1", desc="信息分析流程")
addElement(dom, root, 'p', type="type1", desc=des5)
addElement(dom, root, "pic", name="10X单细胞生信分析流程 图2", type="type1", desc=" ", path="src/images/tu2.png")

# 4 测序数据及其质量评估
des6 = """
高通量测序（如 Illumina HiSeq PE125/PE150）下机得到的原始图像文件经 CASAVA 碱基识别转化为测序读段（Sequenced Reads），
以 FASTQ 格式存储。FASTQ 是一种存储 生物序列及相应质量值的常用文本格式，格式如下。
10X Genomics测序数据每个样本的数据包含I1，R1，R2。I1存储了index信息；R1即read1，
28bp 为细胞 barcode 和 UMI 信息。R2即read2。使用fastqc软件对每个样本的read2数据做质控分析.
"""
addElement(dom, root, 'h1', name="测序数据及其质量评估", type="type1", desc="测序数据及其质量评估")
addElement(dom, root, 'p', type="type1", desc=des6)
addElement(dom, root, "pic", name="FASTQ格式文件示意图 图3", type="type1", desc=" ", path="src/images/tu3.png")

## 4.1 测序碱基质量值
des7 = """
碱基质量值（Quality Score）是碱基识别（Base Calling）出错的概率的整体映射。通常使用的碱基质量值Q公式<a href=\"#ref1\">[3]</a>为：图4。
其中P为碱基识别出错的概率。下表给出了碱基质量值与碱基识别出错的概率的对应关系：
"""
des8 = """
注：碱基质量值越高表明碱基识别越可靠，准确度越高。比如，对于碱基质量值为Q20的碱基识别，100个碱基中有1个会识别出错，以此类推。
"""
des9 = """
注：横坐标为Reads的碱基位置，纵坐标为碱基质量值，蓝线为质量值平均值，箱线的上下线为质量值的极值。
"""
addElement(dom, root, 'h2', name="测序碱基质量值", type="type1", desc="测序碱基质量值")
addElement(dom, root, 'p', type="type1", desc=des7)
addElement(dom, root, "pic", name="碱基质量值Q公式 图4", type="type1", desc=" ", path="src/images/tu4.png")
addElement(dom, root, "pic", name="碱基质量值与碱基识别出错的概率的对应关系 图5", type="type1", desc=des8, path="src/images/tu5.png")
addElement(dom, root, 'p', type="type1", desc='测序数据碱基质量值分布如下图:')
addElement(dom, root, "pic", name="碱基质量分布图 图6", type="type1", desc=des9, path="src/images/tu6.png")

## 4.2 测序碱基含量分布
des10 = """
碱基类型分布检查用于检测有无AT、GC分离现象。理论上，G和C、A和T的含量每个测序循环上应分别相等，
且整个测序过程稳定不变，呈水平线。由于Reads 5’端的前几个碱基为随机引物序列存在一定的偏好性，
因此会在碱基分布图中出现前端波动较大的现象。
"""
addElement(dom, root, 'h2', name="测序碱基含量分布", type="type1", desc="测序碱基含量分布")
addElement(dom, root, 'p', type="type1", desc=des10)
addElement(dom, root, "pic", name="测序碱基含量分布 图7", type="type1", desc=" ", path="src/images/tu7.png")
addElement(dom, root, 'p', type="type1", desc='GC含量分布图如下：')
addElement(dom, root, "pic", name="GC含量分布 图8", type="type1", desc=" ", path="src/images/tu8.png")

## 4.3 原始数据质量评估
des11 = """
#注：SampleID：样本名；Reads：Raw Data中pair-end Reads总数；BaseSum：Raw Data总碱基数；
GC(%):Raw Data GC含量，即Raw Data中G和C两种碱基占总碱基的百分比；Q20(%):
Raw Data质量值大于或等于20的碱基所占的百分比;Q30(%):Raw Data质量值大于或等于30的碱基所占的百分比。
"""
f1 = codecs.open("Rawdata.txt", "wU", "utf-8")
f1.write('\t'.join(['SampleID', 'ReadSum', 'BaseSum', 'GC(%)', 'N(%)', 'Q20(%)', 'Q30(%)' + '\n']))

addElement(dom, root, 'h2', name="原始数据质量评估", type="type1", desc="原始数据质量评估")
addElement(dom, root, 'p', type="type1", desc='该项目各样品数据产出统计见下表：')

for i in range(int(sampleNum)):
    i += 1
    sampleName = 'Sample' + str(i)
    ReadSum = random.randint(310949813, 390949813)
    BaseSum = random.randint(93906843526, 116055237834)
    GC = random.randint(37, 41)
    N = 0
    Q20 = random.randint(81, 83)
    Q30 = random.randint(72, 76)
    f1.write('\t'.join([str(sampleName), str(ReadSum), str(BaseSum), str(GC), str(N), str(Q20), str(Q30) + '\n']))
f1.close()
addElement(dom, root, 'table', name="原始数据产出统计", type="full", desc=des11, path="Rawdata.txt")

# 5 下游分析结果
# 5.1 CellRanger数据统计及定量
des12 = """
注：sampleID：样本ID；
Number of Reads：reads总数；
Valid Barcodes：有效的10X Barcode比例 ；
Sequencing Saturation：测序饱和度；
Q30 Bases in Barcode：Barcode序列中质量值大于或等于30的碱基所占的百分比；
Q30 Bases in RNA Read：reads中质量值大于或等于30的碱基所占的百分比；
Q30 Bases in Sample Index：sample index中质量值大于或等于30的碱基所占的百分比；
Q30 Bases in UMI：UMI序列中质量值大于或等于30的碱基所占的百分比。
"""
addElement(dom, root, 'h2', name="该项目各样品测序数据产出统计见下表：", type="type1", desc="该项目各样品测序数据产出统计见下表：")
f1 = codecs.open("cellranger_results.txt", "wU", "utf-8")
f1.write('\t'.join(['sampleID', 'Number of Reads', 'Valid Barcodes', 'Sequencing Saturation', 'Q30 Bases in Barcode',
                    'Q30 Bases in RNA Read', 'Q30 Bases in Sample Index', 'Q30 Bases in UMI\n']))
for i in range(int(sampleNum)):
    i += 1
    sampleName = 'Sample' + str(i)
    Number_of_Reads = random.randint(310949813, 454288867)
    Valid_Barcodes = str(random.randint(93, 97)) + '%'
    Sequencing_Saturation = str(random.randint(61, 80)) + '%'
    Q30_Bases_in_Barcode = str(random.randint(94, 96)) + '%'
    Q30_Bases_in_RNA_Read = str(random.randint(88, 96)) + '%'
    Q30_Bases_in_Sample_Index = '--'
    Q30_Bases_in_UMI = str(random.randint(92, 96)) + '%'
    f1.write('\t'.join([str(sampleName), str(Number_of_Reads), str(Valid_Barcodes), str(Sequencing_Saturation),
                        str(Q30_Bases_in_Barcode), str(Q30_Bases_in_RNA_Read), str(Q30_Bases_in_Sample_Index),
                        str(Q30_Bases_in_UMI + '\n')]))
f1.close()
addElement(dom, root, 'table', name="CellRanger分析序列统计", type="full", desc=des12, path="cellranger_results.txt")

des13 = """
10x Genomics Chromium 系统利用8通道的微流体 “双十字” 交叉系统，
将含有barcode 的凝胶珠（Gel Beads）、细胞和酶的混合物、油三者混合，
形成 GEMs（油包水的微体系）。理想情况下凝胶珠中的细胞只有1个，
但实验过程中也存在空细胞或者是2个甚至多个细胞的情况。
并且当细胞发生死亡或者裂解时，细胞中将含有大量的线粒体基因。
CellRanger分析时并未对有效细胞进行筛选，因此后续使用Seurat软件做进一步分析。
主要包括细胞过滤、数据标准化处理、PCA分析、t-SNE分析、UMAP分析、细胞聚类、差异Marker Gene分析等。
"""
addElement(dom, root, 'h1', name="下游分析结果", type="type1", desc="下游分析结果")
addElement(dom, root, 'p', type="type1", desc=des13)

## 5.2 数据比对
des14 = """
采用10X Genomics官方软件CellRanger对测序数据进行比对及定量。
CellRanger将read2通过STAR软件比对到参考基因组上。
基于 STAR<a href=\"#ref1\">[4]</a> 的比对结果，结合参考数据集（gtf／gff 文件）里的信息，统计基因组上各个区域的reads覆盖信息，
可以得到比对到外显子，内含子，基因间区的比例信息，做为数据质控的参考指标。
将reads 既比对到已知转录本的外显子上又在同一条链上的作为比对到转录本上的依据。
如果该reads比对到已知的单个基因上，将reads称为唯一比对到转录组上。
只有比对到转录本上的reads才能作为 UMI 计数。
STAR是一款RNA-Seq数据分析常用的分段比对工具，可以用来发现外显子的连接以及融合现象，其基本工作原理主要分成两步：
种子序列的寻找，以及聚类／连接／打分。下图为其原理示意图：
"""
des15 = """
注：最大可比对标签（Maximum Mappable Prefix, MMP）的寻找用于发现（a）外显子连接处，（b）错配，（c）多聚 A 尾， 或者接头，或者低质量尾。
"""
addElement(dom, root, 'h2', name="数据比对", type="type1", desc="数据比对")
addElement(dom, root, 'p', type="type1", desc=des14)
addElement(dom, root, "pic", name="STAR 比对原理 图9", type="type1", desc=des15, path="src/images/tu9.png")
addElement(dom, root, 'p', type="type1", desc='CellRanger分析比对结果统计如下表：')

f1 = codecs.open("mapping.txt", "wU", "utf-8")
f1.write('\t'.join(['sampleID', 'Reads Mapped to Genome', 'Reads Mapped Confidently to Genome',
                    'Reads Mapped Confidently to Intergenic Regions',
                    'Reads Mapped Confidently to Intronic Regions',
                    'Reads Mapped Confidently to Exonic Regions',
                    'Reads Mapped Antisense to Gene',
                    'Reads Mapped Confidently to Transcriptome',
                    'Fraction Reads in Cells\n']))

des16 = """
注：sampleID：样本ID；
Reads Mapped to Genomes：比对到参考基因组上的Reads在总Reads中占的百分比；
Reads Mapped Confidently to Genome：比对到参考基因组并得到转录本GTF信息支持的Reads在总Reads中占的百分比；
Reads Mapped Confidently to Intergenic Regions：比对到基因间区域的Reads在总Reads中占的百分比；
Reads Mapped Confidently to Intronic Regions：比对到内含子区域的Reads在总Reads中占的百分比；
Reads Mapped Confidently to Exonic Regions：比对到外显子区域的Reads在总Reads中占的百分比；
Reads Mapped Antisense to Gene：比对到基因反义链的Reads在总Reads中占的百分比；
Reads Mapped Confidently to Transcriptome：比对到已知参考转录本的Reads在总Reads中占的百分比；
Fraction Reads in Cells：比对到参考基因且来源于高质量细胞的Reads在总Reads中占的百分比；
"""
for i in range(int(sampleNum)):
    i += 1
    sampleName = 'Sample' + str(i)
    a1 = str(random.randint(89, 93)) + '%'
    a2 = str(random.randint(82, 87)) + '%'
    a3 = str(random.randint(2, 3)) + '%'
    a4 = str(random.randint(60, 65)) + '%'
    a5 = str(random.randint(1, 2)) + '%'
    a6 = str(random.randint(1, 2)) + '%'
    a7 = str(random.randint(55, 63)) + '%'
    a8 = str(random.randint(77, 83)) + '%'

    f1.write('\t'.join([sampleName, a1, a2, a3, a4, a5, a6, a6, a8 + '\n']))
f1.close()

addElement(dom, root, 'table', name="CellRanger分析比对结果统计", type="full", desc=des16, path="mapping.txt")

## 5.3 细胞聚类分析
des16 = """
筛选得到有效细胞后，进一步对数据做标准化处理。10X Genomics通过CellRanger软件定量出的表达数据是一个 M*N 的矩阵（行是基因，列是细胞）。
细胞数目通常可达几百上千个。对这样的矩阵聚类计算量极大。因此在对细胞进行聚类之前，
需要先对数据进行降维。目前常用的降维算法包括：PCA(Principal Component Analysis)、
t-SNE(t-Distributed Stochastic Neighbor Embedding)、
UMAP(Uniform Manifold Approximation and Projection)。
PCA降维是一种线性降维方法。运用方差分解，将高维的数据映射到低维的空间中表示。
t-SNE与UMAP是一种用于探索高维数据的非线性降维算法。最终基于SNN聚类算法对细胞进行聚类。
"""
addElement(dom, root, 'h2', name="细胞聚类分析", type="type1", desc="细胞聚类分析")
addElement(dom, root, 'p', type="type1", desc=des16)
addElement(dom, root, "pic", name="Tsne 图10", type="type1", desc=des15, path="src/images/tu10.png")
addElement(dom, root, "pic", name="Umap 图11", type="type1", desc=des15, path="src/images/tu11.png")

## 5.4 差异基因功能富集分析
des16 = """
GO数据库是GO组织（Gene Ontology Consortium）于2000年构建的一个结构化的标准生物学注释系统，
旨在建立基因及其产物知识的标准词汇体系，适用于各个物种。GO注释系统是一个有向无环图，包含三个主要分支，
即：生物学过程（Biological Process），分子功能（Molecular Function）和细胞组分（Cellular Component）。
对每个cluster的差异基因集，采用ClusterProfiler对基因分别进行生物学过程，分子功能和细胞组分的富集分析。
富集分析采用超几何检验方法来寻找与整个基因组背景相比显著富集的GO条目。对富集结果得到的Term采用绘制柱状图气泡图等进行可视化。
"""

des17 = """
在生物体内，不同的基因产物相互协调来行使生物学功能，对差异表达基因的通路（Pathway）注释分析有助于进一步解读基因的功能。
KEGG（Kyoto Encyclopedia of Genes and Genomes）是系统分析基因功能、基因组信息数据库，
它有助于研究者把基因及表达信息作为一个整体网络进行研究。作为是有关Pathway的主要公共数据库(Kanehisa,2008），
KEGG提供的整合代谢途径(pathway)查询，包括碳水化合物、核苷、氨基酸等的代谢及有机物的生物降解，
不仅提供了所有可能的代谢途径，而且对催化各步反应的酶进行了全面的注解，包含有氨基酸序列、PDB库的链接等等，
是进行生物体内代谢分析、代谢网络研究的强有力工具。
对差异表达基因KEGG的注释结果按照KEGG中通路类型进行分类，分类图如下图所示：
"""

addElement(dom, root, 'h2', name="差异基因功能富集分析", type="type1", desc="差异基因功能富集分析")
addElement(dom, root, 'p', type="type1", desc=des16)
addElement(dom, root, "pic", name="GO 图12", type="type1", desc='', path="src/images/tu12.png")
addElement(dom, root, 'p', type="type1", desc=des17)
addElement(dom, root, "pic", name="GO 图13", type="type1", desc='', path="src/images/tu13.png")

## 5.5 细胞拟时序分析
des17 = """
拟时间序列（Pseudotime）简称拟时序，它是一种测量一个细胞在一个生物学过程如细胞分化中进展程度的方法，
它是对一个生物学过程的抽象，描述了一个细胞当前状态和时间轨迹的起始状态之间的最短路径的距离。
拟时分析适用于发育生物学中的发育轨迹研究，或者肿瘤微环境中免疫细胞状态的变化研究等。
Monocle<a href=\"#ref1\">[5]</a>是一款利用反向图形嵌入（reversed graph embedding） 
的机器学习方法，根据细胞基因表达量信息，预测、重构出细胞发育分化轨迹的软件。
"""
addElement(dom, root, 'h2', name="细胞拟时序分析", type="type1", desc="细胞拟时序分析")
addElement(dom, root, 'p', type="type1", desc=des17)
addElement(dom, root, "pic", name="细胞状态(stat)的拟时轨迹图 图14", type="type1", desc='', path="src/images/tu14.png")

## 5.6 细胞亚群的类型鉴定
des18 = """
上述细胞聚类分析是根据细胞之间的相似性将相似度最高的一群细胞识别为一个亚群，但是得到的细胞亚群并没有生物学意义。单细胞注释主要包含两种思路，
一种是将未知细胞类型的细胞数据跟已知类型的细胞数据进行相关性分析，从而对未知的表达数据进行细胞类型鉴定。
另一种方法是根据marker gene进行鉴定，将聚类得到的亚群的marker 基因跟数据库中已知的细胞类型的marker基因进行比较。
SingleR<a href=\"#ref1\">[6]</a>是基于斯皮尔曼相关性，对scRNA-seq数据实现自动化注释的一个软件。百迈客使用 SingleR 对细胞类群进行自动注释（只针对人、小鼠）。
细胞类型鉴定的过程即复杂又繁琐，自动注释只能提供一个参考。
"""
addElement(dom, root, 'h2', name="细胞亚群的类型鉴定", type="type1", desc="细胞亚群的类型鉴定")
addElement(dom, root, 'p', type="type1", desc=des18)
addElement(dom, root, "pic", name="细胞类型分布统计图 图15", type="type1", desc='', path="src/images/tu15.png")
addElement(dom, root, "pic", name="cluster annotations 图16", type="type1", desc='', path="src/images/tu16.png")

# 参考文献：
ref_list = addElement(dom, root, "ref_list", name="参考文献", type="type1", desc="")
addElement(dom, ref_list, "ref", id="1",
           name="Cell RangerTM R Kit Tutorial: Secondary Analysis on 10x GenomicsTM Single Cell 3’ RNA-seq PBMC Data.",
           link="https://cf.10xgenomics.com/supp/cell-exp/cellrangerrkit-PBMC-vignette-knitr-1.1.0.pdf")
addElement(dom, ref_list, "ref", id="2",
           name="Hao*, Hao*, et al., Integrated analysis of multimodal single-cell data bioRxiv 2020.",
           link="https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1")
addElement(dom, ref_list, "ref", id="3",
           name="Ewing B, Hillier L, Wendl MC, Green P. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. [Genome Research Italic] 1998. 8 (3): 175-185.",
           link="http://tour.biocloud.net/article/v1/into/articleDetail/9521921")
addElement(dom, ref_list, "ref", id="4",
           name="Dobin A, Davis C A, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner[J]. Bioinformatics, 2013, 29(1): 15-21.",
           link="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/")
addElement(dom, ref_list, "ref", id="5",  # 15
           name="Xiaojie Qiu, Qi Mao, et al. Reversed graph embedding resolves complex single-cell developmental trajectories.",
           link="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/")
addElement(dom, ref_list, "ref", id="6",  # 16
           name="Aran, Looney, Liu et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology (2019).",
           link="https://www.nature.com/articles/s41590-018-0276-y")

#################################################output my xml object####################################################
domcopy = dom.cloneNode(True)
Indent(domcopy, domcopy.documentElement)
f = open('./' + args.name + '.xml', 'wb')
writer = codecs.lookup('utf-8')[3](f)
domcopy.writexml(writer, encoding='utf-8')
domcopy.unlink()
f.close()
os.system(
    '/share/nas2/genome/biosoft/Python/2.7.8/bin/python  ' + Bin + '/xml2HtmlConverter.py  -i ' + args.name + '.xml' + ' -o ./' + " -l 5 -t " + Bin + '/src ')
