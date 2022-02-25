#!/usr/bin/env perl
my $BEGIN=time();
use strict;
use warnings;
use newPerlBase;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($od,$sample_dir);
$od||="./";
GetOptions(
	"help|?" =>\&USAGE,
	"sample_dir:s"=>\$sample_dir,
) or &USAGE;
&USAGE unless ($sample_dir);
my $Q=33;

my @fq1=glob("$sample_dir/*_S1_L001_R1_001.fastq.gz");
my @fq2=glob("$sample_dir/*_S1_L001_R2_001.fastq.gz");
my @fq;
push @fq,@fq1;
push @fq,@fq2;

foreach my $fq(@fq){
	my $name=basename $fq;
	if($fq=~/fq\.gz/ || $fq=~/fastq\.gz/){
		$name=~s/\.gz$//;
		&runOrDie("gunzip -c $fq > $od/$name");
	}
}

for my $i (0..$#fq1){
	my $f=basename $fq1[$i];
	$f=~s/\_S1_L001_R1_001.fastq.gz//;
	my $a=(glob("$od/${f}_S1_L001_R1_001.fastq"))[0];
	my $b=(glob("$od/${f}_S1_L001_R2_001.fastq"))[0];
	my $cmd =  "fastq_qc_stat -Q $Q -a $a -b $b -f $f  -q 45";
	&runOrDie($cmd);
}

#############################################################
sub USAGE {
	my $usage=<<"USAGE";
Program:
Contact:Liuxs <Liuxs\@biomarker.com.cn>
Description:
Usage:
	Options:
	--sample_dir	sample_dir
	--h	Help
USAGE
	print $usage;
	exit;
}









