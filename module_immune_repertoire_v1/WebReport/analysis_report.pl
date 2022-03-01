#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long::Descriptive;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
print "Program Starts Time:$Time_Start\n";
my $version="1.0";
#######################################################################################
my ($opt, $usage) = describe_options(
        '%c %o --in <infile> --out <outdir> ',
	[ 'cfg1:s', "data_cfg, forced", { required => 1} ],
	[ 'cfg2:s', "detail_cfg, forced", { required => 1} ],
	[ 'sample_filter:s', "dataFilter.singleSample  directory" ],
	[ 'sample_analysis:s', "directory" ],
	[ 'single_enrichment:s', "Anno_enrichment.Annoenrichment" ],
	[ 'single_sample_ppi:s', "single_ppi.sample_ppi" ],
	[ 'single_tf_analysis:s', "single_TF_analysis.TF_analysis" ],
	[ 'single_typeanno:s', "directory" ],
	[ 'integrated_result:s', "data_integrated.integrated_result" ],
	[ 'allcluster_statistic:s', "integrated_analysis.cluster_statistic" ],
	[ 'groupdiff_statistic:s', "group_analysis.groupDiff_statistic" ],
	[ 'group_enrichment:s', "group_enrichment.Annoenrichment" ],
	[ 'cluster_enrichment:s', "cluster_enrichment.Annoenrichment" ],
	[ 'cluster_ppi:s', "cluster_ppi.all_ppi" ],
	[ 'cluster_tf_result:s', "TF_analysis.all_TF_analysis" ],
	[ 'allsample_trace:s', "TraceAnalysis.all_cell_trace" ],
	[ 'cell_typeanno:s', "CellAnno.allAnno" ],
	[ 'cell_cycle:s', "drow_cellcycle.cell_cycle_heatmap" ],
	[ 'od|o:s', "outdir, default ./", { default => "./" } ],
        [ 'help|h', "print help" ],
);

print($usage->text), exit if $opt->help;
#######################################################################
my ($data_cfg,$detail_cfg) = ($opt->cfg1,$opt->cfg2);
my ($sample_filter,$sample_analysis,$single_enrichment,$single_sample_ppi,$single_tf_analysis,$single_typeanno) = ($opt->sample_filter,$opt->sample_analysis,$opt->single_enrichment,$opt->single_sample_ppi,$opt->single_tf_analysis,$opt->single_typeanno);
my ($integrated_result,$allcluster_statistic,$groupdiff_statistic,$group_enrichment,$cluster_enrichment,$cluster_ppi,$cluster_tf_result,$allsample_trace,$cell_typeanno,$cell_cycle) = ($opt->integrated_result,$opt->allcluster_statistic,$opt->groupdiff_statistic,$opt->group_enrichment,$opt->cluster_enrichment,$opt->cluster_ppi,$opt->cluster_tf_result,$opt->allsample_trace,$opt->cell_typeanno,$opt->cell_cycle);
my $od = abs_path($opt->od);

#######################################################################
my $Web_Report = $od;
mkdirOrDie($Web_Report);
my %config;
my %h_sample;
open(CFG1,$data_cfg) or die;
my $num=1;
while(<CFG1>){
        chomp;
        if($_=~/^Sample.*sc$/){
                my $sample=(split /\t|\s+/,$_)[1];
		$h_sample{$sample}="BMK_${num}_$sample";






                $num++;
        }
}
close CFG1;
my %h_group;
$num=1;
open(CFG2,$detail_cfg) or die;
while(<CFG2>){
        chomp;
        next if($_=~/^#/);
        next if($_=~/^\s*$/);
        my ($key,$value)=(split /\t/,$_)[0,1];
        $config{$key}=$value;
        if($_=~/group/){
                my $vs=(split /\t|\s+/,$_)[1];
                $h_group{$vs}="BMK_$num\_$vs";
		$num++;
        }
}
close CFG2;

my $seurat_dir="$Web_Report/BMK_3_seurat_analysis";
	&show_log("Begin BMK_3_seurat_analysis ...");
	mkdirOrDie("$seurat_dir");
	#&readme("$Readme_dir/readme_seurat.pdf", "$seurat_dir");
	mkdirOrDie("$seurat_dir/BMK_1_DataFilter");
	runOrDie("cp $sample_filter/singleSample/*/*.png $seurat_dir/BMK_1_DataFilter || echo skipping!!");
	runOrDie("cp $sample_filter/singleSample/*/*.pdf $seurat_dir/BMK_1_DataFilter || echo skipping!!");

	#单样本的基础分析，降维聚类等
	if(defined $sample_analysis){
		mkdirOrDie("$seurat_dir/BMK_2_Analysis"); # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis
		my @counts_files=glob("$sample_analysis/*/*.All_cell_counts.xls");
		foreach my $counts_file (@counts_files){
			my $filename = basename($counts_file);
			my $dirpath = dirname($counts_file); # ../inputs/analysed/GC03-sc
			my $sample = $1 if($filename=~/(.*).All_cell_counts.xls/);
			my $sample_dir = "$seurat_dir/BMK_2_Analysis/$h_sample{$sample}";
			mkdirOrDie($sample_dir) unless (-d $sample_dir);  # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc
			my $Cluster_dir = "$sample_dir/BMK_1_Clusters";
			mkdirOrDie($Cluster_dir) unless (-d $Cluster_dir);  # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_1_Clusters
			runOrDie("cp $dirpath/*.xls $Cluster_dir || echo skipping!!");
			runOrDie("cp $dirpath/reduction/* $Cluster_dir || echo skipping!!");

			my $diff_analysis = "$sample_dir/BMK_2_Diff_anlaysis";
			mkdirOrDie($diff_analysis) unless (-d $diff_analysis);  # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_2_Diff_anlaysis
			my $clusterFeature = "$diff_analysis/BMK_1_clusterFeature";
			mkdirOrDie($clusterFeature) unless (-d $clusterFeature); # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_2_Diff_anlaysis/BMK_1_clusterFeature
			runOrDie("cp $dirpath/clusterDiff/* $clusterFeature || echo skipping!!");
			runOrDie("cp $dirpath/*.cluster_deg.stat.xls $clusterFeature || echo skipping!!");
			
			my $topmarker = "$diff_analysis/BMK_2_top10_marker";
			mkdirOrDie($topmarker) unless (-d $topmarker); # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_2_Diff_anlaysis/BMK_2_top10_marker
			runOrDie("cp $dirpath/top10/* $topmarker || echo skipping!!");
		}
	}

	# 单样本功能富集分析，以及样本间差异分析
	if(defined $single_enrichment){
		foreach my $single_enrichment_dir (split /,/,$single_enrichment){ # ../inputs/GC03-sc.Anno_enrichment,../inputs/G1Pr-sc.Anno_enrichment
			my $enrich_dirname = (split /\//,$single_enrichment_dir)[-1];
			my $sample=$1 if($enrich_dirname =~ /(.*).Anno_enrichment/);
			my $diff_analysis = "$seurat_dir/BMK_2_Analysis/$h_sample{$sample}/BMK_2_Diff_anlaysis";
			`mkdir -p $diff_analysis` unless (-d $diff_analysis);
			my $Anno_enrichment = "$diff_analysis/BMK_3_Anno_enrichment";
			mkdirOrDie($Anno_enrichment) unless (-d $Anno_enrichment); # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_2_Diff_anlaysis/BMK_3_Anno_enrichment
			runOrDie("cp -r $single_enrichment_dir/* $Anno_enrichment || echo skipping!!!");
		}
	}

	# 单样本 ppi分析，以及样本间差异分析
	if(defined $single_sample_ppi){
		foreach my $single_ppi_dir (split /,/,$single_sample_ppi){ # ../inputs/G1Po-sc.ppi_result,../inputs/G1Pr-sc.ppi_result
			my $ppi_dirname = (split /\//,$single_ppi_dir)[-1];
			my $sample=$1 if($ppi_dirname=~ /(.*).ppi_result/);
			my $diff_analysis = "$seurat_dir/BMK_2_Analysis/$h_sample{$sample}/BMK_2_Diff_anlaysis";
                        `mkdir -p $diff_analysis` unless (-d $diff_analysis);
			my $single_ppi = "$diff_analysis/BMK_4_PPI";
			mkdirOrDie($single_ppi) unless (-d $single_ppi); # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_2_Diff_anlaysis/BMK_4_PPI
			runOrDie("cp -r $single_ppi_dir/* $single_ppi || echo skipping!!!");
		}
	}

	# 单样本 tf分析，以及样本间差异分析
	if(defined $single_tf_analysis){
		foreach my $single_tf_dir (split /,/,$single_tf_analysis){ #../inputs/G1Po-sc.TFBS_Analysis,../inputs/G1Pr-sc.TFBS_Analysis
			my $tf_dirname = (split /\//,$single_tf_dir)[-1];
			my $sample=$1 if($tf_dirname=~ /(.*).TFBS_Analysis/);
			my $diff_analysis = "$seurat_dir/BMK_2_Analysis/$h_sample{$sample}/BMK_2_Diff_anlaysis";
			`mkdir -p $diff_analysis` unless (-d $diff_analysis);
			my $single_tf = "$diff_analysis/BMK_5_TF_analysis/";
			my $single_tf_TFBS = "$single_tf/TFBS_Analysis";
			mkdirOrDie($single_tf) unless (-d $single_tf); # Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_2_Diff_anlaysis/BMK_5_TF_analysis
			mkdirOrDie($single_tf_TFBS) unless (-d $single_tf_TFBS);# Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_2_Diff_anlaysis/BMK_5_TF_analysis/TFBS_Analysis
			runOrDie("cp -r $single_tf_dir/* $single_tf_TFBS || echo skipping!!! || echo skipping!!");
		}
	}

	# 单个样本细胞注释
	if(defined $single_typeanno){
		foreach my $single_typeAnno_dir (split /,/,$single_typeanno){ #../inputs/G1Po-sc.cell_typeAnno,../inputs/G1Pr-sc.cell_typeAnno
			my $anno_dirname = (split /\//,$single_typeAnno_dir)[-1];
			my $sample=$1 if($anno_dirname=~ /(.*).cell_typeAnno/);
			my $sample_dir = "$seurat_dir/BMK_2_Analysis/$h_sample{$sample}";
			`mkdir -p $sample_dir` unless (-d $sample_dir);
			my $s_typeanno = "$sample_dir/BMK_3_cell_typeAnno/";
			mkdirOrDie($s_typeanno) unless (-d $s_typeanno);# Web_Report/BMK_3_seurat_analysis/BMK_2_Analysis/BMK_1_GC03-sc/BMK_3_cell_typeAnno/
			runOrDie("cp -r $single_typeAnno_dir/*cluster_annotation* $s_typeanno || echo skipping!!");
		}
	}

	# 整合的基础结果，包括降维聚类等，可以没有
	unless($integrated_result cmp "None"){
		mkdirOrDie("$seurat_dir/BMK_3_Integrated"); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated
		my $AllData = "$seurat_dir/BMK_3_Integrated/BMK_1_AllData";
		mkdirOrDie($AllData) unless (-d $AllData); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_1_AllData
		runOrDie("cp $integrated_result/base/* $AllData || echo skipping!!");
		my $pca = "$seurat_dir/BMK_3_Integrated/BMK_2_PCA";
		mkdirOrDie($pca) unless (-d $pca); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_2_PCA/
		runOrDie("cp $integrated_result/reduction/ElbowPlot* $pca || echo skipping!!");
		runOrDie("cp $integrated_result/reduction/PCA* $pca || echo skipping!!");
		my $Cluster = "$seurat_dir/BMK_3_Integrated/BMK_3_Cluster";
		mkdirOrDie($Cluster) unless (-d $Cluster); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_3_Cluster
		runOrDie("cp $integrated_result/reduction/tsne* $Cluster || echo skipping!!");
		runOrDie("cp $integrated_result/reduction/umap* $Cluster || echo skipping!!");
		runOrDie("cp $integrated_result/All_ncells_clusters* $Cluster || echo skipping!!");
	}

	# cluster_diff_integrated，整合数据差异分析，可以没有
	unless($allcluster_statistic cmp "None"){
		mkdirOrDie("$seurat_dir/BMK_3_Integrated") unless (-d "$seurat_dir/BMK_3_Integrated");
		my $MarkerGene = "$seurat_dir/BMK_3_Integrated/BMK_4_MarkerGene";
		mkdirOrDie($MarkerGene) unless (-d $MarkerGene); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_4_MarkerGene
		my $Statistics = "$MarkerGene/BMK_1_Statistics";
		mkdirOrDie($Statistics) unless (-d $Statistics);# Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_4_MarkerGene/BMK_1_Statistics
		runOrDie("cp $allcluster_statistic/statistic/* $Statistics || echo skipping!!");
		my $topgene = "$MarkerGene/BMK_2_top10_marker";
		mkdirOrDie($topgene) unless (-d $topgene);# Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_4_MarkerGene/BMK_2_top10_marker
		runOrDie("cp $allcluster_statistic/top10/* $topgene || echo skipping!!");
	}

	# cluster_diff_integrated.Anno_enrichment，整合数据差异分析结果的功能富集分析，可以没有
	unless($cluster_enrichment cmp "None"){
		my $MarkerGene = "$seurat_dir/BMK_3_Integrated/BMK_4_MarkerGene";
		`mkdir -p MarkerGene` unless (-d $MarkerGene);
		my $cluster_enrich = "$MarkerGene/BMK_3_Anno_enrichment";
		mkdirOrDie($cluster_enrich) unless (-d $cluster_enrich);# Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_4_MarkerGene/BMK_3_Anno_enrichment
		runOrDie("cp -r $cluster_enrichment/* $cluster_enrich || echo skipping!!");
	}

	# ../inputs/cluster_diff_integrated.ppi_result，整合数据差异分析结果的ppi分析，可以没有
	unless($cluster_ppi cmp "None"){
		my $MarkerGene = "$seurat_dir/BMK_3_Integrated/BMK_4_MarkerGene";
                `mkdir -p MarkerGene` unless (-d $MarkerGene);
                my $cluster_ppi_dir = "$MarkerGene/BMK_4_PPI";
		mkdirOrDie($cluster_ppi_dir) unless (-d $cluster_ppi_dir);# Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_4_MarkerGene/BMK_4_PPI
		runOrDie("cp -r $cluster_ppi/* $cluster_ppi_dir || echo skipping!!");
	}

	#../inputs/cluster_diff_integrated.TFBS_Analysis，整合数据差异分析结果的tf分析，可以没有
	unless($cluster_tf_result cmp "None"){
		my $MarkerGene = "$seurat_dir/BMK_3_Integrated/BMK_4_MarkerGene";
                `mkdir -p MarkerGene` unless (-d $MarkerGene);
                my $cluster_tf_dir  = "$MarkerGene/BMK_5_TF_analysis";
		mkdirOrDie($cluster_tf_dir) unless (-d $cluster_tf_dir);# Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_4_MarkerGene/BMK_5_TF_analysis
		my $TFBS_dir = "$cluster_tf_dir/TFBS_Analysis";
		mkdirOrDie($TFBS_dir) unless (-d $TFBS_dir); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_4_MarkerGene/BMK_5_TF_analysis/TFBS_Analysis
		runOrDie("cp -r $cluster_tf_result/* $TFBS_dir || echo skipping!!");
	}

	#../inputs/sample_diff_integrated，整合的差异分析结果，可以没有
	unless($groupdiff_statistic cmp "None"){
                my $group_analysis="$seurat_dir/BMK_3_Integrated/BMK_5_Group_Anlysis";
                mkdirOrDie("$group_analysis") unless(-d $group_analysis); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_5_Group_Anlysis
		my @groups_statistic = glob("$groupdiff_statistic/statistic/*.cluster0.all_featuregene.xls");
                foreach my $g_statistic(@groups_statistic){
                        (my $group=basename $g_statistic)=~s/.cluster0.all_featuregene.xls//g;
                        my $group_dir="$group_analysis/$h_group{$group}";
                        mkdirOrDie("$group_dir") unless(-d $group_dir);	# Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_5_Group_Anlysis/BMK_4_CG01-sc_vs_PC01-sc_PC02-sc
                        my $clusterFeature = "$group_dir/BMK_1_Statistics";
                        mkdirOrDie("$clusterFeature") unless(-d $clusterFeature); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_5_Group_Anlysis/BMK_4_CG01-sc_vs_PC01-sc_PC02-sc/BMK_1_Statistics
                        runOrDie("cp -r $groupdiff_statistic/statistic/$group.cluster* $clusterFeature || echo skipping!!");
			my $topmarker = "$group_dir/BMK_2_Topmarker";
			mkdirOrDie("$topmarker") unless(-d $topmarker); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_5_Group_Anlysis/BMK_4_CG01-sc_vs_PC01-sc_PC02-sc/BMK_2_Topmarker
			runOrDie("cp -r $groupdiff_statistic/topmarker/$group.cluster* $topmarker || echo skipping!!");
                }
        }

	# ../inputs/GBJ-sc_vs_HWX-sc.Anno_enrichment，差异分组的功能富集，可以没有
	unless($group_enrichment cmp "None"){
		foreach my $group_enrich (split /,/,$group_enrichment){
			my $groupfile = (split /\//,$group_enrich)[-1];
			my $group = $1 if($groupfile=~/(.*).Anno_enrichment/);
			my $group_analysis="$seurat_dir/BMK_3_Integrated/BMK_5_Group_Anlysis";
        	        mkdirOrDie("$group_analysis") unless(-d $group_analysis); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_5_Group_Anlysis
			my $group_dir= "$group_analysis/$h_group{$group}";
			mkdirOrDie("$group_dir") unless(-d $group_dir); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_5_Group_Anlysis/BMK_4_CG01-sc_vs_PC01-sc_PC02-sc
			my $enrich = "$group_dir/BMK_3_enrichment";
			mkdirOrDie("$enrich") unless(-d $enrich); # Web_Report/BMK_3_seurat_analysis/BMK_3_Integrated/BMK_5_Group_Anlysis/BMK_4_CG01-sc_vs_PC01-sc_PC02-sc/BMK_3_enrichment
			runOrDie("cp -r $group_enrich/* $enrich  || echo skipping!!");
		}
	}

	# 整合数据的注释，可以没有
	unless($cell_typeanno cmp "None"){
		my $cell_type = "$seurat_dir/BMK_4_cell_typeAnno";
		mkdirOrDie("$cell_type") unless(-d $cell_type); # Web_Report/BMK_3_seurat_analysis/BMK_4_cell_typeAnno
		runOrDie("cp -r $cell_typeanno/*cluster_annotation* $cell_type || echo skipping!!");
	}

	#../inputs/analysed_integrated.cell_cycle，整合数据的周期分析，可以没有
	unless($cell_cycle cmp "None"){
		my $cell_cycle_dir = "$seurat_dir/BMK_5_Cell_Cycle";
		mkdirOrDie("$cell_cycle_dir") unless(-d $cell_cycle_dir); # Web_Report/BMK_3_seurat_analysis/BMK_5_Cell_Cycle
		runOrDie("cp -r $cell_cycle/* $cell_cycle_dir || echo skipping!!");
	}

	#../inputs/analysed_integrated.cell_trace，整合数据的轨迹分析，可以没有
	unless($allsample_trace cmp "None"){
		my $allsample_trace_dir = "$seurat_dir/BMK_6_trace_analysis";
		mkdirOrDie("$allsample_trace_dir") unless(-d $allsample_trace_dir); # Web_Report/BMK_3_seurat_analysis/BMK_6_trace_analysis
		runOrDie("cp -r $allsample_trace/*analysed_integrated* $allsample_trace_dir || echo skipping!!!");
	}	
	


#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print STDOUT "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################

sub readme{
	my($readme,$dir)=@_;
	$readme=abs_path($readme);
	runOrDie("cp $readme $dir/readme.pdf");
}


sub T{
        my $id=shift;
        $id = '"'.$id.'"';
        return $id;
}

sub dirmk{
	my $dir=shift;
	&T($dir);
	`rmdir $dir` if(-d $dir);
	mkdirOrDie("$dir");
}

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}
