#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd;
use newPerlBase;

my ($file,$anno_dir,$odir);

GetOptions(
	"help|?" =>\&USAGE,
	"file|f:s" =>\$file,
	"od|o:s" =>\$odir,
	"dir|d:s"=>\$anno_dir,
) or &USAGE;
&USAGE unless($file || $anno_dir);
$odir||="./";
`mkdir $odir` unless(-d $odir);
my %Class=(
        "J" => [1,"Translation, ribosomal structure and biogenesis"],
        "A" => [2,"RNA processing and modification"],
        "K" => [3,"Transcription"],
        "L" => [4,"Replication, recombination and repair"],
        "B" => [5,"Chromatin structure and dynamics"],
        "D" => [6,"Cell cycle control, cell division, chromosome partitioning"],
        "Y" => [7,"Nuclear structure"],
        "V" => [8,"Defense mechanisms"],
        "T" => [9,"Signal transduction mechanisms"],
        "M" => [10,"Cell wall/membrane/envelope biogenesis"],
        "N" => [11,"Cell motility"],
        "Z" => [12,"Cytoskeleton"],
        "W" => [13,"Extracellular structures"],
        "U" => [14,"Intracellular trafficking, secretion, and vesicular transport"],
        "O" => [15,"Posttranslational modification, protein turnover, chaperones"],
        "C" => [16,"Energy production and conversion"],
        "G" => [17,"Carbohydrate transport and metabolism"],
        "E" => [18,"Amino acid transport and metabolism"],
        "F" => [19,"Nucleotide transport and metabolism"],
        "H" => [20,"Coenzyme transport and metabolism"],
        "I" => [21,"Lipid transport and metabolism"],
        "P" => [22,"Inorganic ion transport and metabolism"],
        "Q" => [23,"Secondary metabolites biosynthesis, transport and catabolism"],
        "R" => [24,"General function prediction only"],
        "S" => [25,"Function unknown"],
);


my %Gene;
open IN,"$file" || die;
my $header=<IN>;
chomp $header;
$Gene{'title'}=$header;
while (<IN>) {
        chomp;
        my $id=(split/\s+/,$_)[0];
        $Gene{$id}=$_;
}
close IN;

my %Anno_Gene;
my @Anno_files=glob "$anno_dir/*.anno.txt";
my $ko_anno=glob "$anno_dir/*.ko";
my $pathway_anno=glob "$anno_dir/*.pathway";
my $eggNOG=glob "$anno_dir/*.eggNOG_class.txt";
my $cog=glob "$anno_dir/*.Cog_class.txt";
my $kog=glob "$anno_dir/*.Kog_class.txt";
my $cds_key="Known.longest_transcript.fa";
push @Anno_files,$ko_anno;

my %Anno_Base;

foreach my $file (@Anno_files){
	my $database_key;
        if ($file=~/$cds_key\.Kegg.ko/) {
                $database_key="KEGG";
        }
        elsif ($file=~/$cds_key/ && $file=~/anno/) {
                $file=~m/$cds_key\.(.*).anno.txt/;
                $database_key=$1;
        }
	$Anno_Base{$database_key}=1;
	if ($database_key!~/GO/) {
                open IN,"$file" || die;
                while (<IN>) {
                        chomp;
                        next if (/^$/ || /^\#/);
                        my @annotate=split/\t+/,$_;
                        my $id=$annotate[0];
                        my $anno=$annotate[-1];
                        $Anno_Gene{$id}{$database_key}=$anno;
                }
                close IN;
        }
	else {
                open IN,"$file" || die;
                while (<IN>) {
                        chomp;
                        next if (/^$/ || /^\#/);
                        my ($id,$anno)=split/\t+/,$_,2;
                        $anno=~s/^\d+\t+//;
                        $anno=~s/\t+/\; /g;
                        $Anno_Gene{$id}{$database_key}=$anno;
                }
                close IN;
        }

}

my %Anno_pathway;
if ($pathway_anno) {
                open IN,"$pathway_anno" || die;
                while (<IN>) {
                        chomp;
                        next if (/^$/ || /^\#/);
                        my ($pathway,$pathway_id,$gene_id)=(split/\t+/,$_)[0,1,3];
                        my @gene_id=split/;/,$gene_id;
                        foreach my $gene(@gene_id){
                                push @{$Anno_pathway{$gene}},"$pathway ($pathway_id)";
			}
		}
		close IN;
}

if ($eggNOG) {
    $Anno_Base{eggNOG}=1;

    open NOG,"$eggNOG" || die;
    while (<NOG>) {
        chomp;
        next if (/^$/ || /^\#/) ;
        my @anno=split/\t+/,$_;
	my $classDes="";
        for my $class (split//,$anno[4]){
            next if $class !~/[A-Z]/;
            if(exists $Class{$class}){
                $classDes.=";; $Class{$class}[1]";
            }
            else{
                $classDes.=";; --";
            }
        }
        $classDes=~s/^;; //;
       $Anno_Gene{$anno[0]}{eggNOG}="$anno[4]\t$classDes";
    }
    close NOG;
}

if ($cog) {
    $Anno_Base{COG}=1;

    open COG,"$cog" || die;
    while (<COG>) {
        chomp;
        next if (/^$/ || /^\#/) ;
        my @anno=split/\t+/,$_;
	my $classDes="";
        for my $class (split//,$anno[-3]){
            next if $class !~/[A-Z]/;
            if(exists $Class{$class}){
                $classDes.=";; $Class{$class}[1]";
            }
            else{
                $classDes.=";; --";
            }
        }
        $classDes=~s/^;; //;
       $Anno_Gene{$anno[0]}{COG}="$anno[-3]\t$classDes";
    }
    close COG;
}

if ($kog) {
    $Anno_Base{KOG}=1;

    open KOG,"$kog" || die;
    while (<KOG>) {
        chomp;
        next if (/^$/ || /^\#/) ;
        my @anno=split/\t+/,$_;
	my $classDes="";
        for my $class (split//,$anno[-3]){
            next if $class !~/[A-Z]/;
            if(exists $Class{$class}){
                $classDes.=";; $Class{$class}[1]";
            }
            else{
                $classDes.=";; --";
            }
        }
        $classDes=~s/^;; //;
       $Anno_Gene{$anno[0]}{KOG}="$anno[-3]\t$classDes";
    }
    close KOG;
}

my $name=basename($file);
$name=~s/\.xls//;
open OUT,">$odir/$name.annotation.xls" || die;
print OUT $Gene{'title'};
foreach (sort keys %Anno_Base) {
        if (/eggNOG/) {
                print OUT "\teggNOG_class\teggNOG_class_annotation";
        }
        elsif (/COG/) {
                print OUT "\tCOG_class\tCOG_class_annotation";
        }
    elsif (/KOG/) {
        print OUT "\tKOG_class\tKOG_class_annotation";
    }
        elsif(/KEGG/i){
                print OUT "\tKEGG_annotation\tKEGG_pathway_annotation";
        }
        else {
                my $new_name=$_;
                $new_name=~s/Swissprot/Swiss\-Prot/i;
                $new_name=~s/nr/NR/i;
                $new_name=~s/nt/NT/i;
                $new_name=~s/Cog/COG/i;
                $new_name=~s/Kog/KOG/i;
                print OUT "\t$new_name"."_annotation";
        }
}

print OUT "\n";

foreach my $id (keys %Gene){
	next if($id=~/title/);
	print OUT $Gene{$id};
	foreach (sort keys %Anno_Base) {
	    if (/eggNOG/ && !defined $Anno_Gene{$id}{eggNOG}) {
        	print OUT "\t--\t--";
	    }
	    elsif (/COG/ && !defined $Anno_Gene{$id}{COG}) {
	         print OUT "\t--\t--";
	    }
	    elsif (/KOG/ && !defined $Anno_Gene{$id}{KOG}) {
	        print OUT "\t--\t--";
	    }
	    elsif(/KEGG/i && !defined $Anno_Gene{$id}{KEGG}){
	        print OUT "\t--\t--";
	    }
	    elsif(/KEGG/i && defined $Anno_Gene{$id}{KEGG}){
	        if (defined $Anno_pathway{$id}) {
	            my $pathway=join(";; ",@{$Anno_pathway{$id}});
	            print OUT "\t$Anno_Gene{$id}{$_}\t$pathway";
	        }
	        else{
	            print OUT "\t$Anno_Gene{$id}{$_}\t--";
	        }
		}
	    elsif (defined $Anno_Gene{$id}{$_} and $_!~/KEGG/i) {
	        printf OUT "\t$Anno_Gene{$id}{$_}";
	    }
	    else {
	        printf OUT "\t--";
	    }
	}
	print OUT "\n";
}
close OUT;

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Usage:
  Options:
	-file|f	gene list file
	-dir|d	known gene annotation dir name
	-od|o	outdir name
	-help|h	print usage and exit

USAGE
	print $usage;
	exit;
}





































