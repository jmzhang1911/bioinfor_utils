#!/usr/bin/env perl
my $version="1.0.0";
my $BEGIN=time();
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $programe_dir=basename($0);
my $path=dirname($0);
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;

#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info
# ------------------------------------------------------------------
my ($odir,$HELP,$cfg,$dir);

$odir="./";
GetOptions(
		"cfg:s"=>\$cfg,
		"idir:s"=>\$dir,
		"help"=>\$HELP,
	) or &USAGE;
&USAGE if (defined $HELP);
&USAGE if (!defined $cfg && !defined $dir) ;

my @cluster=glob("$dir/*diff_featuregene.xls");
my $dir_name=(split /\.statistic/,basename($dir))[0];
my %config=();
&readcfg($cfg);
my $ppi=$config{'PPI'};
$odir="$odir/${dir_name}.ppi_result";
`mkdir -p $odir` unless(-d $odir);
for my $l (@cluster) {
	next unless $l=~/[A-Za-z]/;
	my $name = (split /\.diff_featuregene/,basename($l))[0];
	&abstract_interacts($l, "$ppi", "$odir/$name.DEG.detail.txt");
	#&runOrDie("plot_ppi.R --cfg $cfg --indir $dir");
	#&make_sif("$odir/$name.DEG.detail.txt", "$odir/$name.ppi.cytoscapeInput.sif");
	#my $line=`wc -l $odir/$name.ppi.cytoscapeInput.sif|cut -d " " -f 1`;chomp $line;
}
&runOrDie("plot_ppi.R --cfg $cfg --indir $dir");
#`cat $dir/cluster*.gene.xls |grep -v "#" |cut -f 1|sort |uniq >$odir/used_gene.list `;
`cat $dir/*diff_featuregene.xls |grep -v "#" |cut -f 2|sort |uniq >$odir/used_gene.list `;
&abstract_interacts("$odir/used_gene.list", "$ppi", "$odir/ppi_qurey.ppi.detail.txt");
#&make_sif("$odir/ppi_qurey.ppi.detail.txt", "$odir/ppi_qurey.ppi.cytoscapeInput.sif");

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n$programe_dir Start Time :[$Time_Start]\n\n";


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\n$programe_dir End Time :[$Time_End]\n\n";
&Runtime($BEGIN);

#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+
sub abstract_interacts {
    my ($fId, $fIn, $fOut) = @_;

    my %ids;
    open (ID,$fId) or die $!;
    while (<ID>) {
        chomp;
        next if (/^\s*$/ or /^#/);
	my @tmp=split(/\t/, $_);
	my $g;
	if($#tmp<1){
		$g=$tmp[0];
	}else{
		$g=$tmp[1];
	}
        $ids{$g}=1;
    }
    close ID;

    open (IN,$fIn) or die $!;
    open (OUT,">$fOut") or die $!;
    while (<IN>) {
        chomp;
        next if (/^\s*$/);
        if (/^#/) {
            print OUT "$_\n";
        } else {
            my ($k1, $k2) =(split /\t/)[1,4];
            print OUT "$_\n" if ( exists $ids{$k1} && exists $ids{$k2} );
        }
    }
    close IN;
    close OUT;
}

sub make_sif {
    my ($fIn, $fOut) = @_;
    #print "$fIn,$fOut\n";
    my %ppi;
    open (PPI,$fIn) or die $!;
    while (<PPI>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        my ($a,$b,$c) = (split /\t/)[1,2,4];
        ($a,$c) = ($a lt $c) ? ($a,$c) : ($c,$a);
        $ppi{$a}{$c} = $b;
    }
    close PPI;

    open (OUT,">$fOut") or die $!;
    for my $i (sort keys %ppi) {
        for my $j (sort keys %{$ppi{$i}}) {
            print OUT "$i\tpp\t$j\n";
        }
    }
    close OUT;
}

sub readcfg{
    my $cfg=shift;
    open(CFG,$cfg)||die $!;
    while(<CFG>){
        next if($_=~/^#|^\s+$/);
        chomp;my @tmp=split(/\s+|\t/,$_);
        $config{$tmp[0]}=$tmp[1];
    }
    close(CFG);
}

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE 
{
	print <<"	Usage End.";
	Version: $version
	Contact: zhang xuechuan <zhangxc\@biomarker.com.cn>

	Description:
	  -cfg		detail_cfg.cfg
	  -idir		input dir which cluster gene
	  -help
	Usage End.
	exit;
}
