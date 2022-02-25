#!/usr/bin/env perl
#
use strict;
use warnings;
use List::Util qw(max);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#use Algorithm::Combinatorics qw(combinations permutations);
use newPerlBase;
use Cwd 'abs_path';
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
#
my ($od,$cfg,$indir, $id_list);
$od="./";
GetOptions(
	"help|?" =>\&USAGE,
	"indir:s" =>\$indir,
	"id_list:s" =>\$id_list,
	"cfg:s" =>\$cfg,
) or &USAGE;
&USAGE unless ($indir and $cfg);

######
my $dir_name = (split /\.statistic/,basename($indir))[0];
mkdir "$od/${dir_name}.Anno_enrichment" unless (-d "$od/${dir_name}.Anno_enrichment");

my %Number;
my %Symbols;
&GetDEG($indir);
my $column=1;
my %config=();
&readcfg($cfg);
$id_list ||= $config{'Symbol'};
if(defined $id_list){
    &GetList($id_list);
    $column = 2;
}

######
my ($pathway, $go_anno, $go_list, $Known_anno);
$pathway=glob("$config{'Known_anno'}/Known.longest_transcript.fa.Kegg.pathway");
$go_anno=glob("$config{'Known_anno'}/Known.longest_transcript.fa.GO.anno.txt");
$go_list=glob("$config{'Known_anno'}/Known.longest_transcript.fa.GO.list.txt");
$Known_anno=$config{'Known_anno'};

&GetKegg($pathway,"$od/KEGG.info") unless(-e "$od/KEGG.info");
&GetGO($go_anno,"$od/GO.info") unless(-e "$od/GO.info");

open (ANNO,">$od/${dir_name}_enrich.sh") or die $!;
foreach my $cluster (keys %Number){
    my $name=(split /\.diff_featuregene/,basename($cluster))[0];
    if($Number{$cluster} <=10){
		print "the file $cluster DE gene is less than 10.\n";
		`mkdir $od/${dir_name}.Anno_enrichment/$name` unless(-d "$od/${dir_name}.Anno_enrichment/$name");
		`echo "the numbers of $name DE gene is less than 10." >$od/${dir_name}.Anno_enrichment/$name/readme.txt`;
		next;
    }
    if(defined $config{'length'}){
        print ANNO "enrich_analysis.R --deg $cluster --go $od/GO.info --kegg $od/KEGG.info --prefix $name --od $od/${dir_name}.Anno_enrichment/$name --len $config{'length'} --column $column && ";
    }else{
        print ANNO "enrich_analysis.R --deg $cluster --go $od/GO.info --kegg $od/KEGG.info --prefix $name --od $od/${dir_name}.Anno_enrichment/$name --all --column $column && ";
    }
    print ANNO "Reactome_test.R --deg $cluster --prefix $name --od $od/${dir_name}.Anno_enrichment/$name --species $config{'Species'} && ";
    print ANNO "Kegg_map_web.pl -d $cluster -i $Known_anno -k $name -o $od/${dir_name}.Anno_enrichment/$name && ";
    print ANNO "draw_top_GO_graph.pl -i $go_list -deg $cluster -k $name -od $od/${dir_name}.Anno_enrichment/$name && ";
    print ANNO "anno_integrate.pl -f $cluster -d $Known_anno -o $od/${dir_name}.Anno_enrichment/$name && ";
    print ANNO "copy_enrichment_result_dir.pl -i $od/${dir_name}.Anno_enrichment/$name -o $od/${dir_name}.Anno_enrichment/$name \n";
}
close ANNO;

&runOrDie("sh $od/${dir_name}_enrich.sh");

my $END_TIME=time();
my $t=$BEGIN_TIME-$END_TIME;
print "times:$t s";
#######################################################################################
&timeLog("$Script Done");
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
################################################################################################################
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub GetDEG {
	my $dir=shift;
	my @clusters=glob("$dir/*diff_featuregene.xls");
	foreach my $cluster (@clusters){
	    open(T, $cluster) or die $!;
	    my @temp=<T>;
	    close T;
	    $Number{$cluster}=$#temp;
   	}
	my $num=max(values(%Number));
	if($num < 10){
		print "this group's deg gene is less than 10!!!";
		`echo "this group's deg gene is less than 10." > $od/${dir_name}.Anno_enrichment/readme.txt`;
		exit;
	}
}

sub GetList {
	my $list=shift;
	open(LIST,$list) or die $!;
	while(<LIST>){
	    chomp;
            my ($id, $symbol) = (split /\t/, $_)[0..1];
	    if(!defined $Symbols{$id}){
	        $Symbols{$id} = $symbol;
            }
	}
	close LIST;
}

sub GetKegg{
	my ($i,$o)=@_;
	open(IN,$i)||die $!;
	open(OUT,">$o")||die $!;
	while(<IN>){
		chomp;
		next if($_=~/^#/);
		my($path,$ko,$num,$gene,$K)=split(/\t/,$_); 
		my @genes=split(/\;/,$gene);
		foreach my $g(@genes){
			next if($g=~/^$/);
			if(defined $id_list){
			    print OUT "$ko\t$path\t$Symbols{$g}\n";
			}else{
			    print OUT "$ko\t$path\t$g\n";
			}
		}
	}
	close(OUT);
	close(IN);
}

sub GetGO{
	my ($i,$g2go)=@_;
	open(IN,$i)||die $!;
	open(OUT,">$g2go")||die $!;
	while(<IN>){
		chomp;
		next if($_=~/^#/);
		my @tmp=split(/\t/,$_); 
		my $gene=shift @tmp;
		my $num=shift @tmp;
		foreach my $t(@tmp){
			$t=~/(.*?): (.*?) \((GO:.*?)\)/;
			if(defined $id_list){
		   	    print OUT "$3\t$2\t$Symbols{$gene}\t$1\n";
			}else{
			    print OUT "$3\t$2\t$gene\t$1\n";
			}
		}
		print OUT "";
	}
	close(OUT);
	close(IN);
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



sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Usage:
  Options:
  -indir		gene list dir , forced
  -cfg		detail.cfg,contain Ref_seq
  -h		help document

USAGE
	print $usage;
	exit;
}
