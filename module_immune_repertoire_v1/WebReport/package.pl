#!/usr/bin/perl -w
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="2.0.0";
my $indir;
use Cwd 'abs_path';

#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$odir);
#$odir="./";
GetOptions(
                "help|?" =>\&USAGE,
                "in:s"=>\$indir,
		"cfg:s"=>\$cfg,
    	        "od:s" =>\$odir,
                ) or &USAGE;
&USAGE unless ($indir and $odir);
`mkdir $odir` unless(-d $odir);
# ------------------------------------------------------------------
#
# ---indir is program folder---------------------------------------------------------------
$indir = ABSOLUTE_DIR("$indir");
$odir = ABSOLUTE_DIR("$odir");
my @samples;
my @b_samples;
my @t_samples;
open(CFG1,$cfg) or die;
while(<CFG1>){
	chomp;
	if($_=~/^Sample/){
		my ($sample,$type)=(split /\t|\s+/,$_)[1,2];
    if($type eq "sc"){
      push @samples,$sample;
    }
    if($type eq "b"){
      push @b_samples,$sample;
    } 
    if($type eq "t"){
      push @t_samples,$sample;
    }
	}
}

my $back_up="$odir/Backup";
`mkdir $back_up` unless(-d "$back_up");
`mkdir $back_up/single` unless(-d "$back_up/single");
`cp -r $odir/../../inputs/filtered/upload.Rdata  $back_up/single`;
if(-d "$odir/../../inputs/analysed_integrated"){
  foreach my $s(@samples){
    `mkdir -p $back_up/single/$s` unless(-d "$back_up/single/$s");
    `cp -r $odir/../../inputs/filtered/singleSample/$s/single_*.R* $back_up/single/$s`;
  } 
  `mkdir $back_up/integrated` unless(-d "$back_up/integrated");
  `cp -r $odir/../../inputs/analysed_integrated/single_*.R* $back_up/integrated`;
}else{
  foreach my $s(@samples){
    `mkdir -p $back_up/single/$s` unless(-d "$back_up/single/$s");
    `cp -r $odir/../../inputs/analysed/$s/*.single_seruat.R* $back_up/single/$s`;
  } 
}
#`mkdir $back_up/vdj` unless(-d "$back_up/vdj");
#`cp -r $odir/../VDJ_analysis/*_vdj_sc.R* $back_up/vdj`;

system "cd $odir && tar -zcvf Backup.tar.gz Backup";
`rm -rf $odir/Backup`;

my $cellr="$odir/CellRanger";
`mkdir $cellr` unless(-d "$cellr");
foreach my $s(@samples){
  `mkdir -p $cellr/$s` unless(-d "$cellr/$s");
  `cp -r $odir/../../inputs/$s.sc.origin_results/outs/filtered_feature_bc_matrix $cellr/$s`;
  `cp -r $odir/../../inputs/$s.sc.origin_results/outs/raw_feature_bc_matrix $cellr/$s`;
  `cp -r $odir/../../inputs/$s.sc.origin_results/outs/filtered_feature_bc_matrix.h5 $cellr/$s`;
  `cp -r $odir/../../inputs/$s.sc.origin_results/outs/raw_feature_bc_matrix.h5 $cellr/$s`;
  `cp -r $odir/../../inputs/$s.sc.origin_results/outs/cloupe.cloupe $cellr/$s`;
}
foreach my $s(@b_samples){
  `mkdir -p $cellr/$s` unless(-d "$cellr/$s");
  `cp -r $odir/../../inputs/$s.b.origin_results/outs/all_contig_annotations.csv $cellr/$s`;
  `cp -r $odir/../../inputs/$s.b.origin_results/outs/clonotypes.csv $cellr/$s`;
  `cp -r $odir/../../inputs/$s.b.origin_results/outs/consensus_annotations.csv $cellr/$s`;
  `cp -r $odir/../../inputs/$s.b.origin_results/outs/filtered_contig_annotations.csv $cellr/$s`;
  `cp -r $odir/../../inputs/$s.b.origin_results/outs/vloupe.vloupe $cellr/$s`;
}
foreach my $s(@t_samples){
  `mkdir -p $cellr/$s` unless(-d "$cellr/$s");
  `cp -r $odir/../../inputs/$s.t.origin_results/outs/all_contig_annotations.csv $cellr/$s`;
  `cp -r $odir/../../inputs/$s.t.origin_results/outs/clonotypes.csv $cellr/$s`;
  `cp -r $odir/../../inputs/$s.t.origin_results/outs/consensus_annotations.csv $cellr/$s`;
  `cp -r $odir/../../inputs/$s.t.origin_results/outs/filtered_contig_annotations.csv $cellr/$s`;
  `cp -r $odir/../../inputs/$s.t.origin_results/outs/vloupe.vloupe $cellr/$s`;
}

system "cd $odir && tar -zcvf  CellRanger.tar.gz CellRanger";
`rm -rf $cellr`;



my $data_zip = -e "$odir/biomarker_Web_Report.zip";
my $html_zip = -e "$odir/biomarker_htmlReport.zip";

if (! $html_zip) {
    system "cd $indir && zip -q -r $odir/biomarker_htmlReport.zip src/ index.html HTML/";

}

if (! $data_zip) {
    system "cd $indir && zip -q -r $odir/biomarker_Web_Report.zip BMK* HTML index.html src/";

}


#system  "rm -r $indir/BMK_3_mRNA/BMK_7_DEU " if (-d "$indir/BMK_3_mRNA/BMK_7_DEU");
system  "rm -r $indir/template* ";
system "rm -r $indir/Template";
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
    my ($in) = @_;
    my $ret = abs_path($in);
    if (-e $ret ){
        return $ret;
    }else{
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }
}

################################################################################################################

sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   Yang nan <yangn\@biomarker.com.cn> 
Program Date:   2016.03.08
      Modify:   linhj <linhj\@biomarker.com.cn> 
 Modify Date:   2016.03.15

 Description:   This program is used to package the zip file to Need_Data......
       Usage:
        Options:
        -in <dir>   Input directory, this dir must be contain Needed_Data and Web_Report, forced.
        -h          help

USAGE
    print $usage;
    exit;
}
