#!/usr/bin/env perl
use strict;
use Getopt::Long::Descriptive;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
use newPerlBase;

my ($opt,$usage)=describe_options(
	'this program is to analysis TF activity and TFBS %o',
	['all|a=s','all cluster gene\'s expression file,row is geneforced'],
	['cfg|c=s','detail cfg file,forced'],
	['help|h','print usage information and exit']
);


print($usage->text),exit if($opt->help || (!defined $opt->all && !defined $opt->cfg ));
my ($cfg,$all)=($opt->cfg,$opt->all);
my $od="./";
`mkdir $od` unless(-d $od);
my %config;
&readcfg($cfg);

my $name;
if(basename($all) eq "All_cluster_Markergene_avgExp.xls"){
	$name="clusterDiff";
}else{
	$name=(split(/\.All_/,basename($all)))[0];
}

######TFBS
if(exists $config{spe_id}){
	#DEG_TFBS
	$od="$od/${name}.TFBS_Analysis";
	`mkdir -p $od` unless (-d $od);
	&runOrDie("cut -f 1 $all > $od/gene.list");
	`mkdir -p $od/each_DEgeneRes` unless (-d "$od/each_DEgeneRes");
	`mkdir -p $od/DEG_seqLogo` unless (-d "$od/DEG_seqLogo");
	
	my $TFBSdb=$config{'TFBSdb'};
	my @files=();
	my @pdfs=();
    my @pngs=();
	if (-d $TFBSdb) {
		push @files,glob("$TFBSdb/each_geneRes/*predictRes.txt");
	    push @files,glob("$od/each_geneRes/*predictRes.txt");
		push @pdfs,glob("$TFBSdb/seqLogo/*.pdf");
        push @pngs,glob("$TFBSdb/seqLogo/*.png");
        push @pdfs,glob("$od/seqLogo/*.pdf");
        push @pngs,glob("$od/seqLogo/*.png");
	}else{
        push @files,glob("$od/each_geneRes/*predictRes.txt");
		push @pdfs,glob("$od/seqLogo/*.pdf");
        push @pngs,glob("$od/seqLogo/*.png");
	}
	open(LIST,"$od/gene.list")||die $!;
	my %res=();
	my %pdfres=();
    my %pngres=();
	foreach(my $i=0;$i<@files;$i++) {
		my $f=$files[$i];$f=abs_path($f);
		my $gene=basename($f);
       	my $geneID=(split /\_TFBS/,$gene)[0];
		$res{$geneID}=$f;
	}
	foreach(my $i=0;$i<@pdfs;$i++) {
        my $f=$pdfs[$i];$f=abs_path($f);
        my $gene=basename($f);
        my $geneID=(split /\_TFBS/,$gene)[0];
        $pdfres{$geneID}=$f;
    }
	foreach(my $i=0;$i<@pngs;$i++) {
        my $f=$pngs[$i];$f=abs_path($f);
        my $gene=basename($f);
        my $geneID=(split /\_TFBS/,$gene)[0];
        $pngres{$geneID}=$f;
    }
	while(<LIST>){
		chomp;next if(/^#/);
		my $gene=(split /\s+/,$_)[0];
		if(exists $res{$gene}) {
			`cp -r $res{$gene} $od/each_DEgeneRes`;
		}
		if(exists $pdfres{$gene}) {
            `cp -r $pdfres{$gene} $od/DEG_seqLogo`;
        }
        if(exists $pngres{$gene}) {
            `cp -r  $pngres{$gene} $od/DEG_seqLogo`;
        }
	}	
	close(LIST);
	my @DEcheck=glob("$od/each_DEgeneRes/EN*");
    my @DEplots=glob("$od/DEG_seqLogo/EN*.p*");
    if((@DEcheck < 0) || (@DEplots < 0)){
        die "There are something wrong may exits";
    }
}

my @tfbs=glob("$od/each_DEgeneRes/*_TFBS_predictRes.txt");
my $file=shift(@tfbs);
`cat $file  > $od/allGenes_TFBS_predictRes.xls`;
foreach my $tfbs(@tfbs){
	`cat $tfbs | grep -v "^Model_id" >> $od/allGenes_TFBS_predictRes.xls`;
}

my $line=`wc -l $od/allGenes_TFBS_predictRes.xls|cut -b 1`;
`rm $od/allGenes_TFBS_predictRes.xls` if($line==0);
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

