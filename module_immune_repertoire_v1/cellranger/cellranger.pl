#!/usr/bin/env perl
use strict;
use Getopt::Long::Descriptive;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
use newPerlBase;

my ($opt,$usage) = describe_options(
	'mappint and get count data %o',
	['data|d=s','data cfg'],
	['fqdir|f=s','fqstq directory'],
	['detail|t=s','detail cfg'],
	['od|o=s','outdir'],
	['help|h=s','print usage and exit']
);

print($usage->text),exit if(defined $opt->help || !defined $opt->fqdir || !defined $opt->detail);

my($data,$od,$detail)=($opt->data, $opt->od,$opt->detail);
my $fqdir = $opt->fqdir;
$od||="./";
`mkdir -p $od`unless(-d $od);
`mkdir $od/work_sh` unless(-d "$od/work_sh") ;
##########################
my %ref_ge;
open(REF,$detail) or die;
while(<REF>){
  chomp;
  next if($_=~/^$|^\s+$/);
  if($_=~/^VDJ_Ref_Genome/){
    my $vdj_ref=(split /\t|\s+/,$_)[1];
    $ref_ge{"VDJ_Ref_Genome"}=$vdj_ref;
  }
  if($_=~/^SC_Ref_Genome/){
    my $sc_ref=(split /\t|\s+/,$_)[1];
    $ref_ge{"SC_Ref_Genome"}=$sc_ref;
  }
}
print "$ref_ge{'VDJ_Ref_Genome'}\n";
print "$ref_ge{'SC_Ref_Genome'}\n";
########################
my $cmd="";
$od=abs_path($od);
my $fqdirs = basename($fqdir);
my @f = split /\./,$fqdirs;
my $type = $f[1];
my $sample = $f[0];
if($type eq "SC" || $type eq "sc"){
#   `mkdir -p $od/SC`unless(-d "$od/SC");
#   $cmd .="cd $od/SC  && /share/nas1/ranjr/Softwares/cellranger-5.0.0/cellranger count --id=$sample --fastqs=$dir --transcriptome=$ref_ge{'SC_Ref_Genome'} --localcores=16  --localvmem=80 --mempercore=10 \n";
	$cmd .="cd $od && cellranger count --id=$sample --fastqs=$fqdir --transcriptome=$ref_ge{'SC_Ref_Genome'} --localcores=16  --localvmem=80 --mempercore=10 \n";
#   push(@SC,"$od/SC/$sample/outs/metrics_summary.csv");
#	`cp $od/$sample/outs/metrics_summary.csv $od/$sample.$type.metrics_summary.csv`;
}elsif($type eq "B" || $type eq "b"){
#    `mkdir -p $od/B`unless(-d "$od/B");
#    $cmd .="cd $od/B  && /share/nas1/ranjr/Softwares/cellranger-5.0.0/cellranger vdj --chain IG --id=$sample --fastqs=$dir --reference=$ref_ge{'VDJ_Ref_Genome'} --localcores=16 --localvmem=80 --mempercore=10 \n";
    $cmd .="cd $od && cellranger vdj --chain IG --id=$sample --fastqs=$fqdir --reference=$ref_ge{'VDJ_Ref_Genome'} --localcores=16 --localvmem=80 --mempercore=10\n";
#    push(@B,"$od/B/$sample/outs/metrics_summary.csv");
#	`cp $od/$sample/outs/metrics_summary.csv $od/$sample.$type.metrics_summary.csv`;
}elsif($type eq "T" || $type eq "t"){
#	`mkdir -p $od/T`unless(-d "$od/T");
#    $cmd .="cd $od/T  && /share/nas1/ranjr/Softwares/cellranger-5.0.0/cellranger vdj --chain TR --id=$sample --fastqs=$dir --reference=$ref_ge{'VDJ_Ref_Genome'} --localcores=16  --localvmem=80 --mempercore=10 \n";
    $cmd .="cd $od  && cellranger vdj --chain TR --id=$sample --fastqs=$fqdir --reference=$ref_ge{'VDJ_Ref_Genome'} --localcores=16  --localvmem=80 --mempercore=10\n";
#    push(@T,"$od/T/$sample/outs/metrics_summary.csv");
#	`cp $od/$sample/outs/metrics_summary.csv $od/$sample.$type.metrics_summary.csv`;
}

open(SH,">$od/work_sh/cellranger.sh") or die;
print SH "$cmd";
close SH;
&run_or_die("sh $od/work_sh/cellranger.sh");

my $result_dir = "$od/$sample.$type.Result";
`mkdir $result_dir`;

if($type eq "SC" || $type eq "sc"){
	`mkdir $od/$sample.$type.analysis`;
	`cp $od/$sample/outs/metrics_summary.csv $od/$sample.$type.metrics_summary.csv`;
	#############################################

	`cp -r $od/$sample/outs/filtered_feature_bc_matrix $result_dir`;
	`ln $od/$sample/outs/web_summary.html $result_dir/$sample.web_summary.html`;
	`ln $od/$sample/outs/possorted_genome_bam.bam $result_dir/$sample.possorted_genome_bam.bam`;
	`ln $od/$sample/outs/possorted_genome_bam.bam.bai $result_dir/$sample.possorted_genome_bam.bam.bai`;
	`cp -r $od/$sample/outs/analysis/* $od/$sample.$type.analysis`;
	print "Error:there is no filtered_feature_bc_matrix, please check!!!" if(!-d "$result_dir/filtered_feature_bc_matrix");
	print "Error: there is no web_summary.html, please check!!!" if(!-e "$result_dir/$sample.web_summary.html");
} else {
	`cp -r $od/$sample/outs/filtered_contig_annotations.csv $result_dir/$sample.$type.filtered_contig_annotations.csv`;
	`cp -r $od/$sample/outs/web_summary.html $result_dir/$sample.web_summary.html`;
	`cp -r $od/$sample/outs/metrics_summary.csv $od/$sample.$type.metrics_summary.csv`;
}

`mv $od/$sample $od/$sample.$type.origin_results`;


##############################################################################qusb
sub qsub{
        my $shfile= shift;
        my $queue="low.q";
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $queue ";
        &run_or_die($cmd);
        return ;
}
sub run_or_die{
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
        return ;
}
sub show_log{
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
        return ($time) ;
}
sub qsubCheck{
        my $sh = shift;
        my @Check_file = glob "$sh*.qsub/*.Check";
        my @sh_file    = glob "$sh*.qsub/*.sh";
        if ( $#sh_file != $#Check_file ) {
                print "Their Some Error Happend in $sh qsub, Please Check..\n";
                die;
        }else {
                print "$sh qsub is Done!\n";
        }
}







