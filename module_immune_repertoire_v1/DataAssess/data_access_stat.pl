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
my ($od,$sample_stats);
$od||="./";
GetOptions(
	"help|?" =>\&USAGE,
	"sample_stats:s"=>\$sample_stats,
) or &USAGE;
&USAGE unless ($sample_stats);

my @file=split/,/,$sample_stats;
foreach my $f(@file){
	`ln $f $od`;
}

my %mapping;
&summary_qc_stat_info($od);
sub summary_qc_stat_info{
	# get parameters
	my ($od) = @_;

	# get all sample stat file
	my @STAT=glob("$od/*.stat");
	my @Paraf;
	foreach my $stat (@STAT) {
		push @Paraf,$stat unless $stat=~/.+cycleQ\.stat/;
	}
	my %raw_data;

	# open summary stat file
	open (OUT,">$od/AllSample_GC_Q.stat");
	print OUT "sampleID\tReadSum\tBaseSum\tGC(%)\tN(%)\tQ20(%)\tQ30(%)\n";
	foreach my $paraf (@Paraf) {
		# get file name and basename
		my $file=basename($paraf);
		if ($file=~/AllSample.data.stat$/){
                        open(OUT1,">$od/Sample_name.txt")||die"can't open $od/Sample_name.txt\n";
                        open(IN,"$od/$file")||die"can't open $od/$file\n";
                        <IN>;
                        while (<IN>)
                        {
                                chomp;
                                my @A=split/\s+/,$_;
                                if (exists $mapping{$A[0]}) {
                                        $raw_data{$mapping{$A[0]}}=$A[1];
                                        print OUT1 $mapping{$A[0]},"\t",$A[0],"\n";
                                }
                        }
                        close(IN);
                        close(OUT1);
                        next;
                }
	}
	foreach my $paraf (@Paraf) {
                # get file name and basename
                 my $file=basename($paraf);
#		print "file: $file\n";
		next if ($file=~/AllSample_GC_Q.stat$/);	# ignore self
        next if ($file=~/AllSample.data.stat$/);
		$file=~s/\.stat$//;

		# get last line information
		open(IN,"$paraf")||die"can't open $paraf\n";
		<IN>;
		my $now;
		while (<IN>)
		{
			chomp;
			my @A=split/\s+/,$_;
			my $str=join"\t", @A[1..5];
			$now="$file\t$str\t$A[7]\n";
		}
		close(IN);
		print OUT $now;
	}
	close(OUT);

    print "statis_output: $od/AllSample_GC_Q.stat\n";
}

#############################################################
sub USAGE {
	my $usage=<<"USAGE";
Program:
	Contact:Liuxs <Liuxs\@biomarker.com.cn>
	Description:
	Usage:
		Options:
			--sample_stats	sample_stats
			--h	Help
USAGE
print $usage;
exit;
}						
