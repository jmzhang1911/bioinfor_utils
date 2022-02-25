#!/usr/bin/env perl
use strict;
use Getopt::Long::Descriptive;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
use newPerlBase;

my ($summary_file);
my ($opt,$usage) = describe_options(
	'mappint and get count data %o',
	['summary_file|sf=s','summary file'],
	['od|o=s','outdir'],
	['help|h=s','print usage and exit']
);

print($usage->text),exit if(defined $opt->help || !defined $opt->summary_file);
my($summary_file,$od)=($opt->summary_file, $opt->od);

$od||="./";
`mkdir -p $od`unless(-d $od);
$od=abs_path($od);
`mkdir $od/work_sh` unless(-d "$od/work_sh") ;
###################################
my @SC;
my @B;
my @T;
my @summary_file = split /,/,$summary_file;
foreach my $s(@summary_file){
	my $s_tmp = basename($s);
	my @tmp = split /\./,$s_tmp;
	if($tmp[1] eq "SC" || $tmp[1] eq "sc"){
		push @SC,$s;
	}elsif($tmp[1] eq "T" || $tmp[1] eq "t"){
		push @T,$s;
	}elsif($tmp[1] eq "B" || $tmp[1] eq "b"){
		push @B,$s;
	}
}

my $SC=join(",",@SC);
my $T=join(",",@T);
my $B=join(",",@B);
my $stat_cmd;
if($SC){
	#$stat_cmd .="/share/nas2/genome/biosoft/R/3.6.1/bin/Rscript $Bin/basic_info_stat_sc.R -f $SC -o $od/statistic \n";
	$stat_cmd .="Rscript $Bin/basic_info_stat_sc.R -f $SC -o $od/statistic \n";
}
if($T){
	#$stat_cmd .="/share/nas2/genome/biosoft/R/3.6.1/bin/Rscript $Bin/basic_info_stat_t.R -f $T -o $od/statistic \n";
	$stat_cmd .="Rscript $Bin/basic_info_stat_t.R -f $T -o $od/statistic \n";
}
if($B){
	#$stat_cmd .="/share/nas2/genome/biosoft/R/3.6.1/bin/Rscript $Bin/basic_info_stat_b.R -f $B -o $od/statistic \n";
	$stat_cmd .="Rscript $Bin/basic_info_stat_b.R -f $B -o $od/statistic \n";
}

open(ST,">$od/work_sh/statistic.sh") or die;
print ST "$stat_cmd";
close ST;
#&qsub("$od/work_sh/statistic.sh");
#&qsubCheck("$od/work_sh/statistic.sh");
&run_or_die("sh $od/work_sh/statistic.sh");

############################################################qsub
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







