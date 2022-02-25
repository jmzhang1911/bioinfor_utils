#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

my ($Project_name,$Project_key,$Diff,$Thred,$fold,$Species,$Genome_version,$refData,$minexp,$mincell,$minUMI,$minGene,$maxGene,$maxpct,$resolution,$SCT);
GetOptions(
	"help|?"=>\&USAGE,
	"Project_name:s"=>\$Project_name,
	"Project_key:s"=>\$Project_key,
	"Species:s"=>\$Species,
	"Genome_version:s"=>\$Genome_version,
	
	"Diff:s"=>\$Diff,
	"Thred:s"=>\$Thred,
	"fold:s"=>\$fold,

	"refData:s"=>\$refData,
	"minexp:s"=>\$minexp,
	"mincell:s"=>\$mincell,
	"minUMI:s"=>\$minUMI,
	"minGene:s"=>\$minGene,
	"maxGene:s"=>\$maxGene,
	
	"maxpct:s"=>\$maxpct,
	"resolution:s"=>\$resolution,
	"SCT:s"=>\$SCT,
) or &USAGE;
&USAGE unless ($Species);

my $od||="./";
my $degcfg||="deg_cfg.cfg";
my $detailcfg||="detail.cfg";
open(DEG,">$degcfg")||die $!;
open(OUT,">$detailcfg")||die $!;
open(Group,">Group_info.stat")||die $!;

print OUT "Project_name\t$Project_name\n";
print OUT "Project_key\t$Project_key\n";
print OUT "Species\t$Species\n";
print OUT "Genome_version\t$Genome_version\n";

if(defined $refData){
	print OUT "refData\t$refData\n";
	print DEG "refData\t$refData\n";
}

if(defined $Thred){
	my ($thred,$value)=split/:/,$Thred;
	print OUT "$thred\t$value\n";
	print DEG "$thred\t$value\n";
}

if(defined $fold){
	print OUT "fold\t$fold\n";
	print DEG "fold\t$fold\n";
}


if(defined $minexp){
	print OUT "minexp\t$minexp\n";
	print DEG "minexp\t$minexp\n";
}

if(defined $mincell){
	print OUT "mincell\t$mincell\n";
	print DEG "mincell\t$mincell\n";
}

if(defined $minUMI){
	print OUT "minUMI\t$minUMI\n";
	print DEG "minUMI\t$minUMI\n";
}

if(defined $minGene){
	print OUT "minGene\t$minGene\n";
	print DEG "minGene\t$minGene\n";
}

if(defined $maxGene){
	print OUT "maxGene\t$maxGene\n";
	print DEG "maxGene\t$maxGene\n";
}

if(defined $maxpct){
	print OUT "maxpct\t$maxpct\n";
	print DEG "maxpct\t$maxpct\n";
}

if(defined $resolution){
	print OUT "resolution\t$resolution\n";
	print DEG "resolution\t$resolution\n";
}

print OUT "topn\t10\n";
print DEG "topn\t10\n";
##print OUT "length\t100\n";

my $num=1;
my @all_group;
if (defined $Diff) {   ##diff s1_s2_vs_s3_s4,s1_vs_s2
  my @com=split/,/,$Diff;
	push @all_group,@com;
    foreach my $com (@all_group){
        print OUT "group$num\t$com\n";
		print Group "group$num\t$com\n";
		open(G,">group$num.txt")||die $!;
		print G "group\t$com\n";
		close G;
		$num=$num+1;
    }
}
close Group;

if(defined $SCT){
	print OUT "SCT\t$SCT\n";
	print DEG "SCT\t$SCT\n";
}



#my $ref_loc="/share/nas1/ranjr/singleCellProject/genome";
my $Genome_loc="/share/nas2/database/genome";
#my $VDJ_loc="/share/nas1/zhengly/pipeline/10x_vdj/ref/";
my $VDJ_loc="/share/nas2/database/genome/10x_vdj/ref/";
my $ref_loc="/share/nas2/database/genome/scRNA_genome/";


if(defined $Species){
  print OUT "SC_Ref_Genome\t$ref_loc/$Species/$Genome_version/\n" ;
  print OUT "VDJ_Ref_Genome\t$VDJ_loc/${Species}_${Genome_version}/\n";
  print OUT "PPI\t$ref_loc/$Species/$Genome_version/Unigene_Annotation/ppi.txt\n";
  print OUT "Known_unigene\t$ref_loc/$Species/$Genome_version/Unigene_Annotation/Known.longest_transcript.fa\n";
  print OUT "Known_anno\t$ref_loc/$Species/$Genome_version/Unigene_Annotation/Result\n";
  print OUT "Symbol\t$ref_loc/$Species/$Genome_version/genes/id_name.list\n";
}

if($Species=~/Homo_sapiens|Mus_musculus|Rattus_norvegicus/){
	my ($TFDB_dir,$TF,$TFDB);
	$TFDB_dir="/share/nas2/database/AnimalTFDB/v3.0";
	$TFDB="$TFDB_dir/${Species}_TF.txt";
	print OUT "TFDB\t$TFDB\n";
	print OUT "TFBSdb\t$Genome_loc/$Species/$Genome_version/Special/${Species}_TFBS_${Genome_version}\n";
	if($Species=~/Homo_sapiens/){print OUT "spe_id\t9606\n";}
	if($Species=~/Mus_musculus/){print OUT "spe_id\t10090\n";}
	if($Species=~/Rattus_norvegicus/){print OUT "spe_id\t10116\n";}
	print OUT "medical\t$Species\n";
	print DEG "medical\t$Species\n";
}

close OUT;
close DEG;


##############################################
sub USAGE {
	my $usage=<<"USAGE";
	Program:
	Version: $version
	Contact:Liuxs <Liuxs\@biomarker.com.cn> 
	Description:
	Usage:
		Options:
		-Project_name 
		-Project_id
		-Diff
		-Thred
		-fold
		-Species
		-Genome_version
		-refData
		-minexp
		-mincell
		-minUMI
		-minGene
		-maxGene
		-maxpct
		-resolution
		-SCT
		-h	Help
USAGE
	print $usage;
	exit;
}
