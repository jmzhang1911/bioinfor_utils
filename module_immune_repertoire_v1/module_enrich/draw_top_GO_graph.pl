#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Encode;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;

my $BEGIN_TIME=time();
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$deg,$key,$od);
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"deg:s"=>\$deg,
	"k:s"=>\$key,
	"od:s"=>\$od,
) or &USAGE;
&USAGE unless ($fIn and $deg and $od);
$key ||= "output";
mkdir $od unless -d $od;
mkdir "$od/topGO_Graph" unless -d "$od/topGO_Graph";

$fIn = abs_path($fIn);
$deg = abs_path($deg);
$od = abs_path($od);
my $fileType = &terminator_convert($fIn);

my $mid_file_list="$od/GO.list_File_0";
for (my $i=0; ;$i++) {
    if (-f $mid_file_list) {
        my $now=$i+1;
        $mid_file_list="$od/GO.list_File_".$now;
        next;
    }
    last;
}

my $scr = "ReadExcel.r";

if ( $fileType eq 'binary' ) {
    mkdir "$od/annoFiles" unless -d "$od/annoFiles";
print "$scr $fIn $od/annoFiles\n";
    &runOrDie("$scr $fIn $od/annoFiles ");
    my $go_list = "$od/annoFiles/GO.list";
    my @go_list_arr = glob ("$go_list*");
    if (@go_list_arr>1) {
        &runOrDie("cat @go_list_arr > $mid_file_list");
    }elsif (@go_list_arr==0) {
        &submitLog("输入文件异常，请检查！",1,1);
    }else {
        &runOrDie("cp $go_list_arr[0] $mid_file_list");
    }
    &runOrDie("sed -i 's/\\s\\+\$//g' $mid_file_list");

}else {
    my $re = &checkGO($fIn,10);
    if ($re) {
        system("cp $fIn $mid_file_list");
    }else {
        &submitLog("输入文件异常，请检查！",1,1);
#        die "invalid input file, Please check!";
    }
}


if (-f "$mid_file_list") {
    print "draw_GO_DAG_map.pl -All_GO $mid_file_list -DEG_list $deg -od $od/topGO_Graph -key $key \n";
	&runOrDie("draw_GO_DAG_map.pl -All_GO $mid_file_list -DEG_list $deg -od $od/topGO_Graph -key $key");
#	`rm $mid_file_list`;
}
 ` rm -r $od/annoFiles $mid_file_list ` ;


my $context =<<"_README_";
+------------------------------------------------------------------------------+
|                          topGO 分析结果说明文档                              |
+------------------------------------------------------------------------------+

目录结构及文件说明：
********************************************************************************
结果文件/
|-- $key.topGO_BP_gene.xls       #topGO_BP结果,包含GO号对应的具体的基因id。
|-- $key.topGO_BP.pdf            #topGO_BP图(PDF格式)
|-- $key.topGO_BP.png            #topGO_BP图(PNG格式)
|-- $key.topGO_BP.xls            #topGO_BP结果,不包含具体的基因id。
|-- $key.topGO_CC_gene.xls       #topGO_CC结果中GO号对应的具体的基因id。
|-- $key.topGO_CC.pdf            #topGO_CC图(PDF格式)
|-- $key.topGO_CC.png            #topGO_CC图(PNG格式)
|-- $key.topGO_CC.xls            #topGO_CC结果,不包含具体的基因id。
|-- $key.topGO_MF_gene.xls       #topGO_MF结果中GO号对应的具体的基因id。
|-- $key.topGO_MF.pdf            #topGO_MF图(PDF格式)
|-- $key.topGO_MF.png            #topGO_MF图(PNG格式)
|-- $key.topGO_MF.xls            #topGO_MF结果,不包含具体的基因id。
|-- topGO.list                   #基因列表，如果不出现在deg列表中，第二列为1
|-- topGO.map                    #每个基因对应的GO号
`-- readme.txt    #该说明文档

_README_

my $win_context = &T($context);
open (README, ">$od/topGO_Graph/readme.txt") or die $!;
#print README $context;
print README $win_context;
close README;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
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

sub T {
    my $text = shift;
#    return decode( 'gb2312', $text );
	my $line = encode("gb2312", decode("utf8", $text));
	return $line;
}

sub terminator_convert () {
    my $in = shift ;
    my $type = ` $Bin/file -L $in ` ;
    if ($type =~ /with (.*) line terminators/) {
        my $t =$1;
        if ($t =~ /CRLF/) {
            system("$Bin/dos2unix $in");
        }elsif ($t =~ /CR/) {
            system("$Bin/mac2unix $in");
        }
    }elsif ( $type =~ /ASCII text/ ) {
        return "text";
    }else {
        $type = ` $Bin/file -Li $in ` ;
        if ( $type =~ /charset=binary/ ) {
            return "binary";
        }
    }
    return "text";
}

sub checkGO (){
    my ($in,$line) = @_;
    open (IN, "$in") or die "error in [$in], $!";
    while (<IN>) {
        next if ( /^#/ || /^\s*$/) ;
        last if ($. >$line) ;
        chomp;
        my @a = split/\s+/;
        shift @a;
        foreach my $id (@a) {
            if ($id =~ /^GO:\d{7}$/) {
                next;
            }else {
                print "unrecognized format:'$id' in line $. of $in\n";
                return 0;
            }
        }
    }
    close IN;
    return 1;
}
################################################################################################################
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
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
  -i     <file>  input file, All_Database_annotation.xls or GO.list.txt ,forced 
  
  -deg   <file>  deg file,forced 
  
  -k     <str>   keywords of output file, [output] 
  
  -od    <file>  output dir,forced 
  
  -h         Help

USAGE
#	print $usage;
#	exit;
    &submitLog("$usage",1,1);
}
