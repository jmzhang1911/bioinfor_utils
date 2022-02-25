#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;


my ($sample_fq_dir);
GetOptions(
	"help|?"=>\&USAGE,
	"sample_fq_dir:s"=>\$sample_fq_dir
) or &USAGE;
&USAGE unless ($sample_fq_dir);
my $od="./";
my @fq2=glob("$sample_fq_dir/*_S1_L001_R2_001.fastq.gz");
my $cmd;
foreach my $fq(@fq2){
  $cmd .= "fastqc -t 10 $fq -o $od\n";
}
runOrDie("$cmd");


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Contact:Liuxs <Liuxs\@biomarker.com.cn>
Description:
Usage:
  Options:
  --sample_fq_dir  <file>
  --h         Help
USAGE
	print $usage;
	exit;
}
