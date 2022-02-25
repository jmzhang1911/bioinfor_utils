#!/usr/bin/perl -w
use strict;
use File::Basename qw(basename dirname);
use FindBin qw($Script $Bin);
use Data::Dumper;
use Getopt::Long::Descriptive;
use Cwd qw(abs_path);
use newPerlBase;

my ($opt,$usage)=describe_options(
'Description: TCR analysis
Version: Version v1.0
Contact: zhengly <zhengly@biomarker.com.cn> 
Program Date:   2020.12.14 %o',
	['input|i=s','samples cellranger result'],
	['od|o=s','outdir'],
	['cfg|c=s','detail config'],
	['help|h','print usage and exit']
);

print($usage->text),exit if(defined $opt->help);
print($usage->text),exit if(!defined $opt->input and !defined $opt->cfg);
my $in=$opt->input;
my $od=abs_path($opt->od);
my $cfg=abs_path($opt->cfg);
my %config;
&readcfg($cfg);

`mkdir -p $od` unless(-d $od);
`mkdir -p $od/input` unless(-d "$od/input");
`mkdir -p $od/combined` unless(-d "$od/combined");
`mkdir $od/work_sh`;

my @input=split(/,/,$in);
foreach my $data(@input){
   my $name=(split(/\/outs\/filtered_contig_annotations.csv/,$data))[0];
   $name=basename($name);
  `cp -r $data $od/input/$name.csv`;
}

open(SH1,">$od/work_sh/TCR_analysis.sh") or die $!;

my $cmd="/share/nas2/genome/biosoft/R/3.6.1/bin/Rscript $Bin/TCR.R -i $od/input -o $od/TRA -t TRA -s $config{'Species'}\n";

$cmd.="/share/nas2/genome/biosoft/R/3.6.1/bin/Rscript $Bin/TCR.R -i $od/input -o $od/TRB -t TRB -s $config{'Species'}\n";

$cmd.="/share/nas2/genome/biosoft/R/3.6.1/bin/Rscript $Bin/Stat.R -i $od/input -o $od -t TCR\n";

if(-e "$od/../../SC_analysis/step3_integrated/integrated/single_analysis.Rda"){
  $cmd.="export PATH=/share/nas2/genome/biosoft/gcc/5.5.0/bin/:\$PATH && export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/gcc/5.5.0/lib:/share/nas2/genome/biosoft/gcc/5.5.0/lib64:/share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/gsl/2.2.1/lib/:\$LD_LIBRARY_PATH && /share/nas2/genome/biosoft/R/4.0.3/bin/Rscript $Bin/scRepertoire.R -i $od/input -o $od/combined -r $od/../../SC_analysis/step3_integrated//integrated/single_analysis.Rda -t T";
}else{
  my @rda=glob("$od/../../SC_analysis/step2_analysis/*/seurat/single_analysis.Rda");
  my $rda=join("",@rda);
  if(-e "$rda"){
    $cmd.="export PATH=/share/nas2/genome/biosoft/gcc/5.5.0/bin/:\$PATH && export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/gcc/5.5.0/lib:/share/nas2/genome/biosoft/gcc/5.5.0/lib64:/share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/gsl/2.2.1/lib/:\$LD_LIBRARY_PATH && /share/nas2/genome/biosoft/R/4.0.3/bin/Rscript $Bin/scRepertoire.R -i $od/input -o $od/combined -r $rda -t T";
  }
}
print SH1 "$cmd";

close SH1;


&qsub("$od/work_sh/TCR_analysis.sh");
&qsubCheck("$od/work_sh/TCR_analysis.sh");



### sub function
sub writeSH{
	my($cmd,$file) = @_;
	`rm -f $file` if(-e $file);
	open(OUT,">$file")||die $!;
        print OUT "$cmd\n";
        close(OUT);
}

sub readcfg{
        my $cfg = shift;
        open(CFG,$cfg)||die $!;
        while(<CFG>){
                next if($_ =~/^#|^\s+$/);
                chomp;my @tmp = split(/\s+/,$_);
                $config{$tmp[0]} = $tmp[1];
        }
        close(CFG);
}

sub qsub{
        my $shfile= shift;
        my $queue="low.q";
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $queue ";
        &runOrDie($cmd);
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