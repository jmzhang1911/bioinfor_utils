#!/usr/bin/env perl
use strict;
use Getopt::Long::Descriptive;
use File::Basename qw(basename dirname);
use FindBin qw($Script $Bin);

my ($opt,$usage)=describe_options(
	'enrichment result %o',
	['idir|i=s','input dir'],
	['od|o=s','output dir'],
	['help|h','print usage and exit']
);

print($usage->text),exit if(defined $opt->help);
print($usage->text),exit if(!defined $opt->idir and !defined $opt->od);

my $idir=$opt->idir;
my $od=$opt->od;

my $bmk_1="$od/BMK_1_Anno";
my $bmk_2="$od/BMK_2_GO_enrichment";
my $bmk_3="$od/BMK_3_KEGG_enrichment";
my $bmk_4="$od/BMK_4_Reactome_enrichment";
`mkdir $bmk_1` unless(-d $bmk_1);
`mkdir $bmk_2` unless(-d $bmk_2);
`mkdir $bmk_3` unless(-d $bmk_3);
`mkdir $bmk_4` unless(-d $bmk_4);

system "mv $idir/*.annotation.xls $bmk_1";

system "mv $idir/*_Biological_Process_*.png $bmk_2";
system "mv $idir/*_Biological_Process_*.pdf $bmk_2";
system "mv $idir/*_Biological_Process_enrich.list $bmk_2";
system "mv $idir/*_Cellular_Component_*.png $bmk_2";
system "mv $idir/*_Cellular_Component_*.pdf $bmk_2";
system "mv $idir/*_Cellular_Component_enrich.list $bmk_2";
system "mv $idir/*_Molecular_Function_*.png $bmk_2";
system "mv $idir/*_Molecular_Function_*.pdf $bmk_2";
system "mv $idir/*_Molecular_Function_enrich.list $bmk_2";
system "mv $idir/*_go_enrich_barplot.png $bmk_2";
system "mv $idir/*_go_enrich_barplot.pdf $bmk_2";
system "mv $idir/topGO_Graph $bmk_2";

system "mv $idir/*KEGG_pathway_enrich* $bmk_3";
system "mv $idir/pathway/kegg_enrichment/*KEGG.png $bmk_3";
system "mv $idir/pathway/kegg_enrichment/*KEGG.svg $bmk_3";
`mkdir $bmk_3/BMK_1_kegg_map` unless(-d "$bmk_3/BMK_1_kegg_map");
system "cp -r $idir/pathway/kegg_map/* $bmk_3/BMK_1_kegg_map/";

system "mv $idir/*_reactome_enrich.list $bmk_4";
system "mv $idir/*_reactome_enrich*.png $bmk_4";
system "mv $idir/*_reactome_enrich*.pdf $bmk_4";

`rm -r $idir/pathway`;
`rm $bmk_2/topGO_Graph/readme.txt`;
`rm $idir/*cluster*_GO_enrich.list`;
