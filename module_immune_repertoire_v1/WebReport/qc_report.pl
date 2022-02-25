#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin $Script);
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use List::Util qw/max min sum/;
use Cwd qw(abs_path);
use Getopt::Std;
use Config::General;
use IO::File;
use XML::Writer;
use Encode;

my ($input,$output,$project,$cfg,$samples_num);

GetOptions(
	"input:s"	=> \$input,
	"output:s"	=> \$output,
	"project:s"	=> \$project,
	"cfg:s"	=> \$cfg,
	"samples_num:s"	=> \$samples_num,
) || &help;
&help unless ($output and $cfg);

$input = abs_path($input);
mkdir($output,0755) unless -d $output;
$output = abs_path($output);
$cfg = abs_path($cfg);
print $input;
#######
unless (-d "$output/Web_Report/QC_Report/"){
	system "mkdir -p $output/Web_Report/QC_Report/";
}

###更改当前工作路径
chdir "$output/Web_Report/" or warn $!;

#################################
open CFG,$cfg or die $!;
my %cfg;
while(<CFG>){
	chomp;
	next if /^$/ || /^\#/;
	my ($id,$path) = split /\s+/;
	$cfg{$id} = $path;
	if ($id =~ m/报告命名|Project_name/g){
		$project = $path;
	}
}
close CFG;

#################################生成xml文件
my %sampleresult;

my $report_time=&gaintime();
my $writer = XML::Writer->new(OUTPUT => 'self');
$writer->xmlDecl("UTF-8");
$writer->startTag('report');
$writer->emptyTag('report_version','value'=>'v1.0');
$writer->emptyTag('report_name','value'=>"$project");
$writer->emptyTag('report_code','value'=>'XXXX');
$writer->emptyTag('report_user','value'=>'XXXX');
$writer->emptyTag('report_user_addr','value'=>'XXXX');
$writer->emptyTag('report_type','value'=>"otherReport");
$writer->emptyTag('report_time','value'=>$report_time);


#########################################################################

open OUV,">$output/Web_Report/QC_Report/qc_data_1.txt" or die $!;
print OUV "指标名称\t评价标准\t是否达标\n";
open OUM,">$output/Web_Report/QC_Report/qc_data_2.txt" or die $!;
print OUM "指标名称\t评价标准\t是否达标\n";
##QC/OUV
my $UseData = 0;
my $UseGC = 0;
my $UseQ20 = 0;
my $UseQ30 = 0;
##cellranger/OUM
my $Valid_Barcodes = 0;
my $Q30Bases_in_Barcode = 0;
my $Q30Bases_in_RNARead = 0;
my $Sequencing_Saturation = 0;
my $Reads_Mapped_to_Genome = 0; 
my $Reads_Mapped_Confidently_to_Transcriptome = 0;
my $FractionReads_in_Cells = 0;
my $Reads_Mapped_to_Any_vdj = 0;
my $Estimated_Number_of_Cells = 0;
my $Q30Bases_in_UMI = 0;
##########################################################################
##碱基质量
my $odir="$output/Web_Report/BMK_1_rawData";
`cp $odir/*base_quality.png $output/Web_Report/QC_Report/` unless(-e "$output/Web_Report/QC_Report/*base_quality.png");
##碱基含量
#`cp $odir/*acgtn.png $output/Web_Report/QC_Report/` unless(-e "$output/Web_Report/QC_Report/*acgtn.png");
`cp $odir/*gc_content.png $output/Web_Report/QC_Report/` unless(-e "$output/Web_Report/QC_Report/*gc_content.png");
##原始数据产出
`cp $odir/AllSample_GC_Q.xls $output/Web_Report/QC_Report/Sample_GC_Q.stat` unless(-e "$output/Web_Report/QC_Report/Sample_GC_Q.stat");
#my $odir2="$output/Web_Report/BMK_3_seurat_analysis/BMK_1_CellsFilter";
#`cp $odir2/cell_stat_info.xls $output/Web_Report/QC_Report/cell_stat_info.stat` unless(-e "$output/Web_Report/QC_Report/cell_stat_info.stat");
##测序数据产出,比对,细胞信息
my @files = glob("$output/Web_Report/BMK_2_cellranger_analysis/BMK_*_*/BMK_1_summary/*_total_*.xls");
foreach (@files){
	my $filename = basename($_);
	`cp $_ $output/Web_Report/QC_Report/` unless (-e "$output/Web_Report/QC_Report/$_");
}
##########################
open GC,"$output/Web_Report/QC_Report/Sample_GC_Q.stat" or die $!;
while(<GC>){
	chomp;
	my @tmp=split /\t/;
	next if ($tmp[0] eq "sampleID");
if($tmp[2] < 90000000000){
                $UseData ++;  
         }
if($tmp[5] < 80){
		$UseQ20 ++;
	}
 if($tmp[6] < 70){
                $UseQ30 ++;
        }
 if($tmp[3] <=30 and $tmp[4] >=40 ){
                $UseGC ++;
        }
}
close(GC);
#########################
my @seqfiles = glob("$output/Web_Report/QC_Report/*_total_seqence_info_stat.xls");
foreach (@seqfiles){
	my $filename = basename($_);
	my $sc_flag;
	if ($filename eq "sc_total_seqence_info_stat.xls"){
		$sc_flag = "T";
	}else{
		$sc_flag = "F";
	}
	open SEQENCE,$_ or die $!;
	while(<SEQENCE>){
    	    chomp;
	        my @tmp=split /\t/;
    	    print "@tmp\n";
	        next if ($tmp[0] eq "sampleID");
			$tmp[2]=~s/%//;
			$tmp[3]=~s/%//;
			$tmp[-3]=~s/%//;
			$tmp[-2]=~s/%//;
			$tmp[-1]=~s/%//;
			if($tmp[2] < 85){
                $Valid_Barcodes ++;
        	}	
			if(($sc_flag eq "T") && ($tmp[3] < 50)){
    	        $Sequencing_Saturation ++;
        	}
			if($tmp[-3] < 90){
                $Q30Bases_in_Barcode ++;
        	}
			 if($tmp[-2] < 85){
                $Q30Bases_in_RNARead ++;
        	}
			if($tmp[-1] < 85){
				$Q30Bases_in_UMI ++;
			}
	}
	close(SEQENCE);
}
########################
open MAP,"$output/Web_Report/QC_Report/sc_total_mapped_info_stat.xls" or die $!;
while(<MAP>){
        chomp;
        my @tmp=split /\t/;
	 	print "@tmp\n";
        next if ($tmp[0] eq "sampleID");
	$tmp[1]=~s/%//;
	$tmp[7]=~s/%//;
	$tmp[8]=~s/%//;

	if($tmp[1] < 90){        
        $Reads_Mapped_to_Genome ++;
    }  
	 if($tmp[7] < 60){
        $Reads_Mapped_Confidently_to_Transcriptome ++; 
    } 
	 if($tmp[8] < 70 ){
        $FractionReads_in_Cells ++; 
    }
}
close(MAP);
#########################
my @vdjfiles = glob("$output/Web_Report/QC_Report/*_total_vdj_info_stat.xls");
foreach (@vdjfiles){
	open VDJ,$_ or die $!;
	 while(<VDJ>){
		chomp;
		my @tmp=split /\t/;
		print "@tmp\n";
		next if ($tmp[0] eq "sampleID");
		$tmp[1]=~s/%//;
		if($tmp[1] < 50){
			$Reads_Mapped_to_Any_vdj ++;
		}
	}
	close VDJ;
}		
########################
my @cellfiles = glob("$output/Web_Report/QC_Report/*_total_cell_info_stat.xls");
foreach (@cellfiles){
	open CELL ,$_ or die $!;
	while(<CELL>){
		chomp;
		my @tmp=split /\t/;
		print "@tmp\n";
		next if ($tmp[0] eq "sampleID");
		$tmp[1]=~s/,//;
		if($tmp[1] < 10){
			$Estimated_Number_of_Cells ++;
		}
	}
	close CELL;
}


##全符合标准则为yes
if ($UseData == 0){
	print OUV "UseData\t\>\=90\t<font color='#00cc00'>yes</font>\n\n";
}else {
	print OUV "UseData\t\>\=90\t<font color='#ff0000'>no</font>\n";
	}

if ($UseGC == 0){
	print OUV "UseGC\t\>\=30\&\<\=40\t<font color='#00cc00'>yes</font>\n\n";
}else {
	print OUV "UseGC\t\>\=30\&\<\=40\t<font color='#ff0000'>no</font>\n";
	}

if ($UseQ20 == 0){
	print OUV "UseQ20\t\>\=80\t<font color='#00cc00'>yes</font>\n\n";
}else {
	print OUV "UseQ20\t\>\=80\t<font color='#ff0000'>no</font>\n";
	}

if ($UseQ30 == 0){
        print OUV "UseQ30\t\>\=70\t<font color='#00cc00'>yes</font>\n\n";
        }else {
        print OUV "UseQ30\t\>\=70\t<font color='#ff0000'>no</font>\n";
        }
##cellranger
if ($Valid_Barcodes == 0){
        print OUM "Valid_Barcodes\t\>\=95\t<font color='#00cc00'>yes</font>\n\n";
}else {
        print OUM "Valid_Barcodes\t\>\=95\t<font color='#ff0000'>no</font>\n";
}
if ($Q30Bases_in_Barcode == 0){
        print OUM "Q30Bases_in_Barcode\t\>\=85\t<font color='#00cc00'>yes</font>\n\n";
}else {
        print OUM "Q30Bases_in_Barcode\t\>\=85\t<font color='#ff0000'>no</font>\n";
}
if ($Q30Bases_in_RNARead == 0){
        print OUM "Q30Bases_in_RNA_Read\t\>\=85\t<font color='#00cc00'>yes</font>\n\n";
}else {
        print OUM "Q30Bases_in_RNA_Read\t\>\=85\t<font color='#ff0000'>no</font>\n";
}
if ($Q30Bases_in_UMI == 0){
		print OUM "Q30Bases_in_UMI\t\>\=85\t<font color='#00cc00'>yes</font>\n\n";
}else {
		print OUM "Q30Bases_in_UMI\t\>\=85\t<font color='#ff0000'>no</font>\n";
}
if ($Estimated_Number_of_Cells == 0){
		print OUM "Estimated_Number_of_Cells\t\>\=10\t<font color='#00cc00'>yes</font>\n\n";
}else {
		print OUM "Estimated_Number_of_Cells\t\>\=10\t<font color='#ff0000'>no</font>\n";
}

if ($Sequencing_Saturation == 0){
        print OUM "Sequencing_Saturation\t\>\=50\t<font color='#00cc00'>yes</font>\n\n";
}else {
        print OUM "Sequencing_Saturation\t\>\=50\t<font color='#ff0000'>no</font>\n";
}
if ($Reads_Mapped_to_Genome == 0){
        print OUM "Reads_Mapped_to_Genome\t\>\=90\t<font color='#00cc00'>yes</font>\n\n";
}else {
        print OUM "Reads_Mapped_to_Genome\t\>\=90\t<font color='#ff0000'>no</font>\n";
}
if ($Reads_Mapped_to_Any_vdj == 0){
		print OUM "Reads_Mapped_to_Any_V\(D\)J_gene\t\>\=50\t<font color='#00cc00'>yes</font>\n\n";
}else {
		print OUM "Reads_Mapped_to_Any_V\(D\)J_gene\t\>\=50\t<font color='#ff0000'>no</font>\n";
}
if ($Reads_Mapped_Confidently_to_Transcriptome == 0){
        print OUM "Reads_Mapped_Confidently_to_Transcriptome\t\>\=60\t<font color='#00cc00'>yes</font>\n\n";
}else {
        print OUM "Reads_Mapped_Confidently_to_Transcriptome\t\>\=60\t<font color='#ff0000'>no</font>\n";
}
if ($FractionReads_in_Cells == 0){
        print OUM "FractionReads_in_Cells\t\>\=70\t<font color='#00cc00'>yes</font>\n\n";
}else {
        print OUM "FractionReads_in_Cells\t\>\=70\t<font color='#ff0000'>no</font>\n";
}
	close OUM;
	close OUV;
##项目概况
$writer->emptyTag('report_abstract','value'=>"<p class='p-abstract'>本页面展示项目结果质控（QC）信息。第一部分是所有质控指标列表，第二部分是每个指标的详细信息说明和结果展示。</p><p class='p-abstract'>项目信息概况：</p><p class='p-abstract'>本项目共分析了$samples_num 个样品。");
$writer->emptyTag('h1','desc'=>'质控指标列表','type'=>"type1",'name'=>'所有质控指标列表');
my $pid++; 
$writer->emptyTag('p','desc'=>"测序数据及其质量评估，包括碱基质量值、碱基含量发布、原始数据质量评估。具体统计结果如下：",'type'=>"type1");
my $tabletype=&judgetable("$output/Web_Report/QC_Report/qc_data_1.txt");
$writer->emptyTag('table','desc'=>'','name'=>"测序数据质量指标","type"=>"$tabletype",'path'=>"/QC_Report/qc_data_1.txt");

$writer->emptyTag('p','desc'=>"cellranger数据统计及定量结果如下；",'type'=>"type1");
$writer->emptyTag('table','desc'=>"",,'name'=>"比对/定量指标","type"=>"full",'path'=>"/QC_Report/qc_data_2.txt");

#$writer->emptyTag('p','desc'=>"Seurat数据统计结果如下；",'type'=>"type1");
#$writer->emptyTag('table','desc'=>"",,'name'=>"细胞质量指标","type"=>"full",'path'=>"$output/Web_Report/QC_Report/qc_data_3.txt");

###详细说明
$writer->emptyTag('h1','desc'=>'每个指标的详细说明和结果展示','type'=>"type1",'name'=>'每个指标的详细说明和结果展示');
##碱基质量 BMK_1_rawData/
$writer->emptyTag('p','desc'=>" Clean data，请根据项目情况，判定数据产出是否达标，其信息统计见下表：",'type'=>"type1");
&piclist("碱基质量分布图","注","$output/Web_Report/QC_Report/*base_quality.png","/QC_Report/");
##碱基含量 BMK_1_rawData/
#$writer->emptyTag('p','desc'=>" Clean data，请根据项目情况，判定数据产出是否达标，其信息统计见下表：",'type'=>"type1");
#&piclist("碱基含量分布图","注","$output/Web_Report/QC_Report/*acgtn.png","/QC_Report");
##测序数据产出统计 
$writer->emptyTag('p','desc'=>" Clean data，请根据项目情况，判定数据产出是否达标，其信息统计见下表：",'type'=>"type1");
&piclist("GC含量分布图","注","$output/Web_Report/QC_Report/*gc_content.png","/QC_Report");
##原始数据产出统计 BMK_1_rawData/
$writer->emptyTag('p','desc'=>"Clean data，请根据项目情况，判定数据产出是否达标，其信息统计见下表：",'type'=>"type1");
$writer->emptyTag('table','desc'=>"注：SampleID：样本名；Reads：Raw Data中pair-end Reads总数；BaseSum：Raw Data总碱基数；GC(%):Raw Data GC含量，即Raw Data中G和C两种碱基占总碱基的百分比；Q20(%):Raw Data质量值大于或等于20的碱基所占的百分比;Q30(%):Raw Data质量值大于或等于30的碱基所占的百分比。",'name'=>"原始数据产出统计结果","type"=>"full",'path'=>"/QC_Report/Sample_GC_Q.stat");
##CellRanger测序数据产出结果
$writer -> emptyTag('p','desc' => '该项目各样品测序数据产出统计见下表：','type'=> "type1");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Number of Reads：reads总数；Valid Barcodes：有效的10X Barcode比例 ；Sequencing Saturation：测序饱和度；Q30 Bases in Barcode：Barcode序列中质量值大于或等于30的碱基所占的百分比；Q30 Bases in RNA Read：reads中质量值大于或等于30的碱基所占的百分比；Q30 Bases in Sample Index：sample index中质量值大于或等于30的碱基所占的百分比；Q30 Bases in UMI：UMI序列中质量值大于或等于30的碱基所占的百分比。",'name'=>"CellRanger分析序列统计","type"=>"full",'path'=>"/QC_Report/sc_total_seqence_info_stat.xls");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Number of Reads：reads总数；Valid Barcodes：有效的10X Barcode比例 ；Q30 Bases in Barcode：Barcode序列中质量值大于或等于30的碱基所占的百分比；Q30 Bases in RNA Read：reads中质量值大于或等于30的碱基所占的百分比；Q30 Bases in Sample Index：sample index中质量值大于或等于30的碱基所占的百分比；Q30 Bases in UMI：UMI序列中质量值大于或等于30的碱基所占的百分比。",'name'=>"B 细胞序列统计","type"=>"full",'path'=>"/QC_Report/b_total_seqence_info_stat.xls");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Number of Reads：reads总数；Valid Barcodes：有效的10X Barcode比例 ；Q30 Bases in Barcode：Barcode序列中质量值大于或等于30的碱基所占的百分比；Q30 Bases in RNA Read：reads中质量值大于或等于30的碱基所占的百分比；Q30 Bases in Sample Index：sample index中质量值大于或等于30的碱基所占的百分比；Q30 Bases in UMI：UMI序列中质量值大于或等于30的碱基所占的百分比。",'name'=>"T 细胞序列统计","type"=>"full",'path'=>"/QC_Report/t_total_seqence_info_stat.xls");

##CellRanger分析比对结果
$writer -> emptyTag('p','desc' => '采用10X Genomics官方软件CellRanger对测序数据进行比对及定量。CellRanger将read2通过STAR软件比对到参考基因组上。基于 S统计基因组上各个区域的reads覆盖信息，可以得到比对到外显子，内含子，基因间区的比例信息，做为数据质控的参考指标。','type'=> "type1");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Reads Mapped to Genomes：比对到参考基因组上的Reads在总Reads中占的百分比；Reads Mapped Confidently to Genome：比对到参考基因组并得到转录本GTF信息支持的Reads在总Reads中占的百分比；Reads Mapped Confidently to Intergenic Regions：比对到基因间区域的Reads在总Reads中占的百分比；Reads Mapped Confidently to Intronic Regions：比对到内含子区域的Reads在总Reads中占的百分比；Reads Mapped Confidently to Exonic Regions：比对到外显子区域的Reads在总Reads中占的百分比；Reads Mapped Antisense to Gene：比对到基因反义链的Reads在总Reads中占的百分比；Reads Mapped Confidently to Transcriptome：比对到已知参考转录本的Reads在总Reads中占的百分比；Fraction Reads in Cells：比对到参考基因且来源于高质量细胞的Reads在总Reads中占的百分比；",'name'=>"CellRanger分析比对结果统计","type"=>"full",'path'=>"/QC_Report/sc_total_mapped_info_stat.xls");

##CellRanger B/T 细胞 V(D)J基因的Reads富集和组装
$writer -> emptyTag('p','desc' => 'CellRanger 在对Reads 富集前，使用 cutadapt去除Read-pairs上的接头和引物序列。随后将 Read-pairs比对到 V(D)J 基因片段上，将比对上的 read 用于后续的组装。组装的过程中，每个 barcode 是独立进行分析的。对过滤后的 read 按照 barcode 分组，其中每个barcode至多有100k的reads用于组装，且仅使用多于 10 个 reads 的 UMI 的 read 用于组装，最后获得 contigs 序列。','type'=> "type1");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Reads Mapped to Any V(D)J Gene：比对到V(D)J基因片段的Reads在总Reads中占的百分比；Reads Mapped to IGH/IGK/IGL：比对到IGH/IGK/IGL基因片段的Reads在总Reads中占的百分比；Median IGH/IGK/IGL UMIs per Cell：每个细胞 IGH/IGK/IGL UMIs中位数；Cells With Productive V-J Spanning Pair：拥有成对的有效 V 基因至 J 基因区的细胞占比；Cells With IGH/IGK/IGL Contig：含有 IGH/IGK/IGL contig的细胞比例；；Cells With CDR3-annotated IGH/IGK/IGL Contig：含有CDR3序列的 IGH/IGK/IGL contig的细胞比例；Cells With V-J Spanning IGH/IGK/IGL Contig：含有跨越V基因的5’端和J基因的3’端 IGH/IGK/IGL contig的细胞比例。",'name'=>"CellRanger B细胞分析统计","type"=>"full",'path'=>"/QC_Report/b_total_vdj_info_stat.xls");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Reads Mapped to Any V(D)J Gene：比对到V(D)J基因片段的Reads在总Reads中占的百分比；Reads Mapped to IGH/IGK/IGL：比对到IGH/IGK/IGL基因片段的Reads在总Reads中占的百分比；Median IGH/IGK/IGL UMIs per Cell：每个细胞 IGH/IGK/IGL UMIs中位数；Cells With Productive V-J Spanning Pair：拥有成对的有效 V 基因至 J 基因区的细胞占比；Cells With IGH/IGK/IGL Contig：含有 IGH/IGK/IGL contig的细胞比例；；Cells With CDR3-annotated IGH/IGK/IGL Contig：含有CDR3序列的 IGH/IGK/IGL contig的细胞比例；Cells With V-J Spanning IGH/IGK/IGL Contig：含有跨越V基因的5’端和J基因的3’端 IGH/IGK/IGL contig的细胞比例。",'name'=>"CellRanger T细胞分析统计","type"=>"full",'path'=>"/QC_Report/t_total_vdj_info_stat.xls");

##CellRanger分析细胞结果
$writer -> emptyTag('p','desc' => '10X scRNA-Seq的基因表达定量，主要基于UMI计数来实现的。通过UMI可以区分一条read是否属于生物学重复还是技术重复，能够有效地去除PCR效应。CellRanger对每个Barcode下的基因去除重复的UMI，统计unique UMI数目即表示细胞基因的表达量。CellRanger分析细胞统计如下表：','type'=> "type1");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Estimated Number of Cells：检测到的细胞数目；Mean Reads per Cell：每个细胞平均Reads数目 ；Median UMI Counts per Cell：每个细胞的UMI中位数；Median Genes per Cell：每个细胞中基因数目的中位数；Total Genes Detected：所有细胞的基因总数。",'name'=>"CellRanger分析细胞结果统计","type"=>"full",'path'=>"/QC_Report/sc_total_cell_info_stat.xls");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Estimated Number of Cells：检测到的细胞数目；Mean Reads per Cell：每个细胞平均Reads数目 ；Median UMI Counts per Cell：每个细胞的UMI中位数；Median Genes per Cell：每个细胞中基因数目的中位数；Total Genes Detected：所有细胞的基因总数。",'name'=>"CellRanger分析细胞结果统计","type"=>"full",'path'=>"/QC_Report/b_total_cell_info_stat.xls");
$writer->emptyTag('table','desc'=>"注：sampleID：样本ID；Estimated Number of Cells：检测到的细胞数目；Mean Reads per Cell：每个细胞平均Reads数目 ；Median UMI Counts per Cell：每个细胞的UMI中位数；Median Genes per Cell：每个细胞中基因数目的中位数；Total Genes Detected：所有细胞的基因总数。",'name'=>"CellRanger分析细胞结果统计","type"=>"full",'path'=>"/QC_Report/t_total_cell_info_stat.xls");

#$writer->endTag('ref_list');
$writer->endTag('report');
open OUT,">:utf8", "$output/Web_Report/configtest_qc.xml";
my $xmlstr=&decorate($writer->to_string);
$xmlstr = Encode::decode("utf-8", $xmlstr);
print OUT $xmlstr;
close(OUT);
$writer->end();
exit(0);



####################
sub judgetable{
	my $file=shift;
	my $line=`wc -l $file`;	chop($line);
	$line=(split(/\s+/,$line))[0];
	my $flag= $line>5 ? "full" : "short";
	return $flag;
}

sub piclist{
	my ($name,$desc,$pics,$path,$num)=@_;	##图名称、图注释、图路径、图分割路径、最大图片数（默认20）
		my @images=glob("$pics");
	return if(scalar(@images)==0);
	$num||=23;
	$pid++;	my $i=0;
	$writer->startTag('pic_list','name'=>"图$pid. $name",'desc'=>"$desc",'type'=>"type1");
	foreach my $s(@images){
		if($i==$num){
			$writer->endTag('pic_list');
			return;
		}
		my $base=basename $s;
		my $dir=dirname $s;
		my $tmp=(split(/$path/,$dir))[1];
		my $new =$path.$tmp;
		$writer->emptyTag('pic','desc'=>"",'name'=>"$base",'type'=>"type1",'path'=>"$new/$base");
		$i++;		
	}
	$writer->endTag('pic_list');
	return;
}
sub help()
{
	my $usage =<< "	EOF";
	
	Usage:
	    -input	input directory
	    -output	output directory
	    -project	project name
	    -cfg	detail cfg
		-samples_num	sample num
	    
	EOF
	print $usage;
	exit;
}

sub decorate{
        my $xmlstr=shift;
        $xmlstr=~s/(<[^>]*)\>/$1\>\n/mgo;
        return $xmlstr;
}
sub gaintime{
        my $timestamp=time;
        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($timestamp);
        my $y = $year + 1900;
        my $m = $mon + 1;
        $timestamp=printf("%d-%02d-%02d",$y,$m,$mday);
        #return $timestamp
}

sub digitize
{
    my $v = shift or return '0';
    $v =~ s/(?<=^\d)(?=(\d\d\d)+$)     #处理不含小数点的情况
            |
            (?<=^\d\d)(?=(\d\d\d)+$)   #处理不含小数点的情况
            |
            (?<=\d)(?=(\d\d\d)+\.)    #处理整数部分
            |
            (?<=\.\d\d\d)(?!$)        #处理小数点后面第一次千分位
            |
            (?<=\G\d\d\d)(?!\.|$)     #处理小数点后第一个千分位以后的内容，或者不含小数点的情况
            /,/gx;
    return $v;
}

sub run_or_die()
{
        my ($cmd) = @_ ;
        my $flag = system($cmd) ;
        if ($flag != 0){
                print "Error: command fail: $cmd\n";
                exit(1);
        }
        return ;
}

sub dpic{
        my ($writer,$desc,$plist)=@_;
        my ($sampid,$i,$bn);
        foreach my $path(@$plist)
        {       if(defined $desc)
                {
                        $i=$desc;
                        $bn=basename($path);
                        $bn=~/(X\d+)\./;
                        $sampid=$1;
                        $i=~s/xxx/$sampid/i;
                }
                $writer->emptyTag('pic','name'=>"",type=>"img-width-normal" ,'desc'=>"",'path'=>"$path");
        }
}

sub dpic2{
        my($writer,$desc,$plist) = @_;
        for     (my $var = 0;$var < @$desc;$var++){
                $writer->emptyTag('pic','name'=>"@$desc[$var]",type => "img-width-normal",'desc'=>"@$desc[$var]",'path'=>"@$plist[$var]");
        }
}
