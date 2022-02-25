#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($All_GO,$DEG_list,$key,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"All_GO:s"=>\$All_GO,
				"DEG_list:s"=>\$DEG_list,
				"key:s"=>\$key,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($All_GO and $DEG_list and $key and $od);

mkdir $od unless -d $od;

my %H;
open (IN,"$DEG_list") or die $!;
while (<IN>) {
	next if /^\#/;
	my @A=split/\s+/,$_;
	#my $gene=$A[0];
	my $gene=$A[0];
	my $val=(@A>=4)?$A[2]:0.01;
	# my $val=(@A>=4)?$A[-3]:0.01;
	#my $val=0.01;
	$H{$gene}=$val;
}
close IN;

open (IN,"$All_GO") or die $!;
open (OUT,">$od/topGO.map") or die $!;
open (OUT2,">$od/topGO.list") or die $!;
print OUT2 "#ID\tValue\n";
my $f=0;
while (<IN>) {
	chomp;
	my ($name,$info)=(split/\t/,$_,2)[0,1];
	next unless defined $info;
	$info=~s/\t/,/g;
	print OUT "$name\t$info\n";
    if (exists $H{$name}) {
        print OUT2 "$name\t$H{$name}\n" ;
        $f =1;
    }else {
        print OUT2 "$name\t1\n";
    }
}
close IN;
close OUT;
close OUT2;

if ($f==0) {
    &submitLog("基因列表和背景基因集不匹配，请检查！",1,1);
}

print "topGO.R -o $od -k $key";
&runOrDie("topGO.R -o $od -k $key");

chdir $od;
&runOrDie("svg2xxx -dpi 300 $key.topGO_BP.svg");
&runOrDie("svg2xxx -dpi 300 $key.topGO_MF.svg");
&runOrDie("svg2xxx -dpi 300 $key.topGO_CC.svg");

` rm topGO.list topGO.map `;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub T {
    my $text = shift;
#    return decode( 'gb2312', $text );
	my $line = encode("gb2312", decode("utf8", $text));
	return $line;
}

sub LogErrorFile {
	my $err_info = shift;
    $err_info= &T($err_info);
	open ERROUT,">$od/topGO.error";
	print ERROUT "$err_info";
	close ERROUT;
	#print "/share/nas1/limh/docker_sge/RunningMessageProducer type code \"$err_info\"\n";
	` /usr/local/bin/RunningMessageProducer 1 1 "$err_info"  `;
}

sub submitLog (){
    my $log =shift;
    my $type = shift;  # 1 for stderr;  5 for stdout.
    my $die = shift; 
    # /usr/local/bin/RunningMessageProducer <type> <code> <desc>
    ` /usr/local/bin/RunningMessageProducer $type 200 "$log"  `;
    if ($die && $die==1) {
        die "$log";
    }
}



################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Usage:
  Options:
  -All_GO    <file>  input file,forced 
  
  -DEG_list  <file>  input file,forced 
  
  -key       <str>   index of output files,forced 
  
  -od        <dir>   output dir,forced 
  
  -h         Help

USAGE
#	print $usage;
#	exit;
    &submitLog("$usage",1,1);
}
