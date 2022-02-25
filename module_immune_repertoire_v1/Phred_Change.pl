#!/usr/bin/env perl
use strict;
use Getopt::Long::Descriptive;
use File::Basename qw(basename dirname);
use FindBin qw($Script $Bin);
use Cwd qw(abs_path);
use newPerlBase;

my ($r1_fqs, $r2_fqs, $sample, $type);
my($opt,$usage) = describe_options(
	'Data preprocessing',
	['r1_fqs|r1=s','r1 fqst1'],
	['r2_fqs|r2=s','r2 fastq2'],
	['sample|s=s','sample name'],
	['type|t=s','type'],
	['od|o=s','outdir name'],	
	['help|h','print usage and exit']
);
print($usage->text),exit if(defined $opt->help || !defined $opt->sample || !defined $opt->r1_fqs || !defined $opt->r2_fqs || !defined $opt->type);

#my $data=abs_path($opt->cfg1);
#my $detail=abs_path($opt->cfg2);
#my $od=abs_path($opt->od);

my $sample=$opt->sample;
my $r1_fqs=$opt->r1_fqs;
my $type=$opt->type;
my $r2_fqs=$opt->r2_fqs;

my $od ||= "./";
`mkdir $od` unless(-d $od);
`mkdir $od/work_sh` unless(-d "$od/work_sh");
#########################################
open(CONF,">$od/data.cfg")||die $!;
my (%data_cfg);
my %All_Q;
my @g1_fqs=split/,/,$r1_fqs; ##fq1
#print "r1,$r1_fqs\n";
my @g2_fqs=split/,/,$r2_fqs; ##fq2
my @r1_fqs;my @r2_fqs;
my @sample=split/,/,$sample; ##sample name
my @type=split /,/,$type; ##data type
##########################################
open OUT,">$od/samples_name.xls" or die $!;
print OUT "#ID\ttype\tread1_name\tread2_name\n";

print CONF "Qphred\t33\n";
for(my $i=0;$i<=$#g1_fqs;$i++){
	#$data_cfg1{$sample[$i]}{$type[$i]}=$g1_fqs[$i];
	#$data_cfg2{$sample[$i]}{$type[$i]}=$g2_fqs[$i];
	print CONF "Sample\t$sample[$i]\t$type[$i]\n";
	my $Fastqs1 = basename($g1_fqs[$i]);
	my $Fastqs2 = basename($g2_fqs[$i]);
	print CONF "fq1\t$sample[$i].$type[$i].fastq/$sample[$i]\_S1_L001_R1_001.fastq.gz\n";
	print CONF "fq2\t$sample[$i].$type[$i].fastq/$sample[$i]\_S1_L001_R2_001.fastq.gz\n";
	print OUT "$sample[$i]\t$type[$i]\t$Fastqs1\t$Fastqs2\n";
	
	`mkdir $od/$sample[$i].$type[$i].fastq` unless(-d "$od/$sample[$i].$type[$i].fastq");
	`ln $g1_fqs[$i] $od/$sample[$i].$type[$i].fastq/$sample[$i]\_S1_L001_R1_001.fastq.gz`;
	#`cp $g1_fqs[$i] $od/$sample[$i].$type[$i].fastq/`;
	#`mv $od/$sample[$i].$type[$i].fastq/$Fastqs1 $od/$sample[$i].$type[$i].fastq/$sample[$i]\_S1_L001_R1_001.fastq.gz`;
	`ln $g2_fqs[$i] $od/$sample[$i].$type[$i].fastq/$sample[$i]\_S1_L001_R2_001.fastq.gz`;
	#`cp $g2_fqs[$i] $od/$sample[$i].$type[$i].fastq/`;
	#`mv $od/$sample[$i].$type[$i].fastq/$Fastqs2 $od/$sample[$i].$type[$i].fastq/$sample[$i]\_S1_L001_R2_001.fastq.gz`;
}
close(CONF);
close(OUT);

######################################### sub function
sub get_fq{
    my $fq=shift @_;
    print "the fq file is $fq\n";
    my $name=$fq;
    if ($fq=~/fq\.gz/ || $fq=~/fastq\.gz/){
        $name=$fq;
    }else{
        die "please check the input fq file !";
    }
    return $name;
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











