#!/usr/bin/perl -w
use strict;
use Getopt::Long::Descriptive;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;

my ($opt,$usage)=describe_options(
'Description: get template file 
Version: v1.0
Contact: Ranjr <ranjr@biomarker.com.cn> 
Program Date: 2020.9.17 %o',
	['idir|i=s','input dir, eg:Web_Report'],
	['od|o=s','out dir,eg:Web_Report'],
	['cfg|c=s','detail config file'],
	['type|t=s','type of report'],
	['help|h','print usage and exit']
);

print($usage->text),exit if(defined $opt->help);
print($usage->text),exit if(!defined $opt->idir && !defined $opt->od);

my $indir=abs_path($opt->idir);
my $odir=abs_path($opt->od);
my $cfg=abs_path($opt->cfg);
my $type=$opt->type;
my %config;
&readcfg($cfg,\%config);

`mkdir $odir` unless (-d $odir);
`rm -rf $odir/Template ` if -d "$odir/Template";
mkdir "$odir/Template" unless -d "$odir/Template";
system "cp $Bin/template/* $odir/Template/";
if($type eq 'm'){
	`mv $odir/Template/file_template_muli.txt $odir/Template/file_template.txt`;
}
if($type eq 's'){
	`mv $odir/Template/file_template_sig.txt $odir/Template/file_template.txt`;
}
my @ftemplates;
open(TEM,"<$odir/Template/file_template.txt") or die $!;
while (<TEM>) {
    chomp;
    next if (/^#/ || /^$/);
    my ($file,$template)=split /\s+/,$_,2;
    my $dir;
        if ($file=~/\|/){
                print "$file\n";
                my ($t_file,$demo)=split /\|/,$file,2;
                my ($t1,$t2)=split /\|/,$template,2;
                my $t=(glob"$indir/$t_file")[0];
                if (-e $t){
                        $dir=$t;
                        $template=$t1;
                }else{
                        $dir="$odir/$demo";
                        $template=$t2;
                }
        }else{$dir=(glob "$indir/$file")[0];}
    my $flag;
    if (defined $dir) {
        if (!-f $dir && !-d $dir ) {
            print "$file is not use in the report \n";
        }elsif (-e $dir){
	    push @ftemplates,"$odir/Template/$template";
	}
    }
}

`cat @ftemplates >$odir/Template/s1_template.txt`;

`cat $odir/Template/s1_template.txt > $odir/template`;

my $cells_stat=glob("$odir/BMK_2_cellranger_analysis/BMK_1_SC/BMK_1_summary/sc_total_cell_info_stat.xls");
chomp(my $sample_num=`wc -l $cells_stat | awk '{print \$1}'`);
$sample_num=$sample_num-1;
my $cells_sum=0;
open(IN,$cells_stat) or die $!;
<IN>;
while(<IN>){
	chomp;
	my $cell=(split /\t/,$_)[1];
	$cell=~s/,|;//;
	$cells_sum=$cells_sum+$cell;
}
close(IN);
my $t_cell=0;
if(-e "$odir/BMK_2_cellranger_analysis/BMK_2_T/BMK_1_summary/t_total_cell_info_stat.xls"){
open(TC,"$odir/BMK_2_cellranger_analysis/BMK_2_T/BMK_1_summary/t_total_cell_info_stat.xls") or die $!;
<TC>;
while(<TC>){
	chomp;
	my $cell=(split /\t/,$_)[1];
	$cell=~s/,|;//;
	$t_cell=$t_cell+$cell;
}
close(TC);
}
my $b_cell=0;
if(-e "$odir/BMK_2_cellranger_analysis/BMK_2_B/BMK_1_summary/b_total_cell_info_stat.xls"){
open(BC,"$odir/BMK_2_cellranger_analysis/BMK_2_B/BMK_1_summary/b_total_cell_info_stat.xls") or die $!;
<BC>;
while(<BC>){
	chomp;
	my $cell=(split /\t/,$_)[1];
	$cell=~s/,|;//;
	$b_cell=$b_cell+$cell;
}
close(BC);
}
my $AllSample_GC_Q=glob("$odir/BMK_1_rawData/AllSample_GC_Q.xls");
my $base_num=0;
my $base_t=0;
my $base_b=0;
open(A,$AllSample_GC_Q) or die $!;
<A>;
while(<A>){
	chomp;
	my ($na,$base)=(split /\t/,$_)[0,2];
	$base=~s/,|;//;
	if($na=~/-sc$/){
		$base_num=$base_num+$base;
	}
	if($na=~/-t$/){
		$base_t=$base_t+$base;	
	}
	if($na=~/-b$/){
		$base_b=$base_b+$base;
	}
}
close A;
my $mean_base=($base_num/$sample_num)/1000;
$mean_base=sprintf "%.2f",$mean_base;
my $mean_b=($base_b/$sample_num)/1000;
$mean_b=sprintf "%.2f",$mean_b;
my $mean_t=($base_t/$sample_num)/1000;
$mean_t=sprintf "%.2f",$mean_t;
my %data_among_words;
$data_among_words{'\$FC'}=$config{'fold'};
$data_among_words{'\$FDR'}=$config{'FDR'} if(defined $config{'FDR'});
$data_among_words{'\$FDR'}=$config{'pvalue'} if(defined $config{'pvalue'});
$data_among_words{'\$sample_num'}=$sample_num;
$cells_sum=&format_figure($cells_sum);
#print("$cells_sum\n");
$data_among_words{'\$cells'}=$cells_sum;
$data_among_words{'\$data'}=$mean_base;
$t_cell=&format_figure($t_cell);
#print("$t_cell\n");
$data_among_words{'\$t_cells'}=$t_cell;
$data_among_words{'\$t_data'}=$mean_t;
$b_cell=&format_figure($b_cell);
#print("$b_cell\n");
$data_among_words{'\$b_cells'}=$b_cell;
$data_among_words{'\$b_data'}=$mean_b;
#exit;
my %ref;
&REF_HASH();
&GET_REF("$odir/Template/s1_template.txt","$odir/Template/reference.txt","$odir/Template/s2_template");

if (-f "$odir/Template/reference.txt"){
    `cat $odir/Template/s2_template $odir/Template/reference.txt > $odir/template`;
}

#`rm $odir/Template/template*.txt`;
#`rm $odir/Template/*.txt`;

###############################################################
sub REF_HASH{
    open(REF,"$odir/Template/template.reference.txt") or die $!;
    <REF>;<REF>;
    while (<REF>) {
        chomp;
        my @tem=split /\t/,$_;
        my $num=scalar @tem;
        if ($num!=3) {
            print "$.\t$num\n";
            warn "must have 3 colunm !!";
        }
        $ref{$tem[0]}{'ref'}=$tem[1];
        $ref{$tem[0]}{'link'}=$tem[2];
        
    }
    close REF;
}
sub GET_REF{
	my ($tem_1,$ref,$tem_2)=@_;
    open(OUT,">$tem_2") or die $!;
    open(OUT1,">$ref") or die $!;
    print OUT1 "》\n参考文献\n";
    open(IN,"$tem_1") or die $!;
    my %count;
    my $ref_num=1;
    while (<IN>) {
        chomp;
        next if (/^#/ || /^$/);
        my @te=$_=~/(\$REF_\d+)/g;
        my $n=scalar @te;
        if ($n==0){
            $_=&Trans($_);
            print OUT "$_\n" ;
        }else{
            foreach my  $key (@te){
                $key=~/\$REF_(\d+)/;
                my $id=$1;
                if (!exists $count{$id}) {
                    $_ =~ s/\$REF_$id/[$ref_num]/;
                    $count{$id}=$ref_num;
                    print OUT1 "$ref_num\t$ref{$id}{'ref'}\t$ref{$id}{'link'}\n";
                    $ref_num++;   
                }elsif (exists $count{$id}){
                    $_ =~ s/\$REF_$id/[$count{$id}]/;
                }
            }
            $_=&Trans($_);
            print OUT "$_\n";
        }    
    }
    close IN;
    close OUT;
    close OUT1;
}

sub Trans{
	my $text = shift;
	my (@data_among_word) = $text =~ /(\$[a-z,A-Z,_,0-9]+)/g;
	if (@data_among_word !=0) {
		for (my $j=0;$j<@data_among_word;$j++) {
			my $data_id = $data_among_word[$j];
			$data_among_word[$j] =~ s/\$/\\\$/;
			$text =~ s/$data_among_word[$j]/$data_among_words{$data_among_word[$j]}/;
		}
	}
	my (@ref_among_word) = $text =~ /\[(\d*)\]/g;

	return $text;
}
sub href {
	my $text =shift;
	my $name = '&lt;a href=&quot;#ref'.$text.'&quot;&gt;['.$text.']&lt;/a&gt;';
	return $name;
}

sub readcfg{
	my ( $cfg_file, $config ) = @_;
	open( CFG, $cfg_file ) or die "$!: $cfg_file\n";
	while (<CFG>) {
        	chomp;
	        s/^\s+//;s/\s+$//;s/\r$//;
	        next if ( /^\s+/ or /^#/ || /^$/ );
	        my ( $key, $value ) = ( split /\s+/ )[0,1];
        	$config->{$key} = $value;
    	}
    	close CFG;
    &log_current_time("detail config done.");
}

sub date_time_format {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );
}

sub log_current_time {
	my ($info) = @_;
	my $curr_time = date_time_format( localtime( time() ) );
	print "[$curr_time] $info\n";
}
sub Integer_Three_Digit{#
	my $interger = shift;
	$interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $interger;
}
sub format_figure{#
	my $figure = shift;
	if (!defined $figure) {
		die;
	}
	if ($figure=~/\./) {
		if ($figure == 100) {
			$figure = 100;
		} else {
		$figure = sprintf("%.2f",$figure);
		}
	}else{
		$figure = Integer_Three_Digit($figure);
	}
	return $figure;
}
