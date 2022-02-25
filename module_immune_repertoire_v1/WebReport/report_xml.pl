#!/usr/bin/env perl
use autodie qw(:all);
use strict;
use warnings;
use Getopt::Long::Descriptive;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd qw (abs_path getcwd);
use File::Spec;
use newPerlBase;
use File::Basename qw (basename dirname);
my $BEGIN_TIME = time ();
my $version = "2.0.0";
use RTF::Writer;
use Encode qw (decode);

###############################################################

my ($opt,$usage)=describe_options(
	'convert xml file %o',
	['indir|i=s','input dir'],
	['cfg|c=s','config file'],
	['temp|t=s','template file'],
	['xml|x=s','the xml file'],
	['bio|b','local xml'],
	['help|h','print usage and exit']
);

print($usage->text),exit if(defined $opt->help);
print($usage->text),exit unless($opt->indir and $opt->cfg and $opt->xml and $opt->temp);

my $dir=$opt->indir;
my $config=$opt->cfg;
my $ftemplate=$opt->temp;
my $xml=$opt->xml;
my $bio=$opt->bio;

my $type="type1";
my $desc_anno="";
$desc_anno=&T($desc_anno);
my $table_type="full";
$table_type=&T($table_type);
my $dir_abs = abs_path ($dir);
`rm -r $dir/HTML` if (-d "$dir/HTML");
mkdir "$dir/HTML" unless -d "$dir/HTML";
my $HTML="$dir_abs/HTML";
my $dir_template = "$dir_abs/Template";

my %data;
my %table_info;
my %pic_info;
my %data_among_words;
my %Ref_info;
my %config;
my %file_all_info; 
my %reference;
my %second_links;

open (IN,"$config") or die $!;
while (<IN>) {
	chomp;
	s/\s+$//;s/\r$//;
	next if (/^#/ || /^$/);
	my ($k,$v) = (split /\s+/,$_)[0,1];
	if ($k =~ /^Project_name/){
		$config{$k}=$v;
	}
    if ($k =~ /^Project_id/){
		$config{$k}=$v;
	}
    if ($k=~/^Ref_seq/){
		$config{$k}=$v;
	}
}
close IN;

my $prefix = $config{Project_name};
chomp (my $user = `whoami`);
my $user_addr = $user."\@biomarker.com.cn";
my $report_version = &T ("v1.0");
my $report_name = &substitute($prefix);
$report_name=&T($prefix);
$report_name = &Trans($report_name);
my $report_code = $config{Project_id};
$report_code=&T($report_code);
$report_code = &Trans($report_code);
my $report_user = &T($user);
my $report_user_addr = &T($user_addr);
my $report_time = &GetDate;
my $alltime =$report_time;
$alltime=&T($alltime);
$alltime=&Trans($alltime);
$report_time = &T($report_time);

# -------------------------------------------------------------------------
# output xml
# -------------------------------------------------------------------------
open (OUT,">$dir_template/xml_1.xml") or die $!;
print OUT "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
print OUT "<report>\n";
print OUT "\t<report_version value=$report_version \/>\n";
print OUT "\t<report_name value=$report_name \/>\n";
print OUT "\t<report_code value=$report_code \/>\n";
print OUT "\t<report_user value=$report_user \/>\n";
print OUT "\t<report_user_addr value=$report_user_addr \/>\n";
print OUT "\t<report_time value=$alltime \/>\n";
print OUT "\t<report_type value=\"lncrna\" \/>\n";
print OUT "\t<report_abstract value='' \/>\n";
open (TEM, $ftemplate) or die $!;
$/="》";
my $i=0;
my $pic_nummber=1;
my $tab_nummber=1;
while (<TEM>){
	chomp;
	next if (/^$/ || /^\#/);
	my (undef ,$format,@context)=split /\n/,$_;
	$i++;
	my $num=@context;
	$format =~ s/\s+$//;
	if ($format =~ /^图片/ ){print OUT &picture_write($format,$type,\@context);}
	elsif ($format =~ /^表格/){print OUT &table_write($format,$type,\@context); }
	elsif ($format =~ /^统计表|附表/){print OUT &stat_table_write($format,$type,\@context); }
	elsif ($format =~ /^HTML/){print OUT &HTML_write($format,$type,\@context); }
	elsif ($format =~ /^示意图/){print OUT &picture_write($format,"img-width-normal",\@context); }
	elsif ($format =~ /^示意表/){print OUT &schematic_table_write($format,$type,\@context); }
	elsif($format =~ /^正文/){print OUT &Text_mul($format,$type,\@context);}
	elsif($format =~ /^参考文献/){print OUT &Reference($format,$type,\@context);}
	elsif ($format =~ /^摘要/) { print OUT &abstract($format,$type,\@context);}
	elsif ($format =~ /^公式/) {print OUT &Formula($format,$type,\@context);}
	elsif ($format =~ /级标题/) {print OUT &Head($format,$type,\@context);	}
	elsif($format=~/文件集合/) {print OUT &file_all_write($format,\@context);}
	elsif($format=~/网页报告/) {print OUT &web_all_write($format,\@context);}
	elsif($format=~/结果目录说明/) {print OUT &Web_Report_directory($format,$type,\@context);}
	
}
print OUT "<\/report>\n";
close TEM;
close OUT;
$/="\n";    ###
system "perl $Bin/xml_add_number.pl --xmlin $dir_template/xml_1.xml --xmlout $dir_template/xml_2.xml ";
if (defined $bio){
	open( XML, "$dir_template/xml_2.xml" ) or die $!;
	open( XML1, ">$xml") or die $!;
	while (<XML>) {
		chomp;
		$_ =~ s/$dir_abs\///g;
		my $xml_temp1 = $_;
		print XML1 "$xml_temp1\n";
	}
	close XML;
	close XML1;
}

system "rm $dir_template/xml_*.xml";
`find $HTML/ -name "*.html" |xargs sed -i '/table-responsive/ s/&lt;br&gt;//g'`;


sub picture_write {
	my ($name,$type,$context) = @_;
	my ($idir,$desc,$anno);
	$desc=$$context[0];
	$desc=~ s/图(\d+)\s/图$pic_nummber / if ($desc!~/公式/);
	$pic_nummber++ if ($desc!~/公式/);
	$anno=(@$context==2)?$$context[1]:" ";
	$idir=(split /\s+/,$name)[1];
	my $file="$dir_abs/$idir";
	my @picts = glob($file);
	my $pict_num = @picts;
	my $content;
	$type = &T($type);
	$anno=&T($anno);
	if (1<$pict_num) {
	
			$desc = &T($desc);
			$content = "\t<pic_list name=$desc type=$type desc=$anno >\n";
			my $i=1;
			foreach my $pict(@picts) {
				if($desc=~/差异表达基因的KEGG通路注释图/ && $i>10){
					last;
				}
				my $pict_name = &T (basename($pict));
				my $path = &T($pict);
				$content = $content."\t\t<pic name=$pict_name desc=$desc_anno path=$path \/>\n";
				$i++;
			}
		$content = $content."\t".'</pic_list>'."\n";
		
	}
	if($pict_num==1) {
		my $path=$picts[0];
		if (-f "$path"){
			my $pict_name = &T($desc);
			$path = &T($path);
			$content="\t<pic name=$pict_name type=$type desc=$anno path=$path \/>\n";	
		}
	}
	return $content;
}
sub table_write {
	my ($name,$type,$context) = @_;
	my ($idir,$desc,$anno);
	$desc=$$context[0];
	my $desc1=$desc;
	$desc1=~s/表\d+//;
	$desc=~s/表(\d+)\s/表$tab_nummber / if ($desc!~/附表/) ;
	$tab_nummber++ if ($desc!~/附表/) ;
	$anno=(@$context==2)?$$context[1]:" ";
	$idir=(split /\s+/,$name)[1];
	my $file="$dir_abs/$idir";
	my @tables = glob ($file);
	my $content;
	$type = "xls";
	$type = &T($type);
	$anno=&T($anno);
	my $table_num = @tables;
	my $txt=$tables[0];
	my $txt_line=`less -S $txt|grep -v "#" |wc -l`; chomp $txt_line;
	next if ($table_num ==0);
	$txt=$tables[1] unless ($txt_line !=1);
	$desc = &T($desc);
	$desc1 = &T($desc1);
	if ($txt_line > 1){
		my $table_name=basename($txt);
		my $txt_tmp=shift;
		if(!-f "$dir/Template/$table_name"){
			#`head -n 6 $txt >$dir/Template/$table_name` if($txt_line > 6);
			#`cat $txt > $dir/Template/$table_name` if($txt_line <= 6);
			`head -n 37 $txt >$dir/Template/$table_name` if($txt_line > 37);
			`cat $txt > $dir/Template/$table_name` if($txt_line <= 37);
			$txt_tmp="$dir/Template/$table_name";
		}else{
			#`head -n 6 $txt >$dir/Template/$table_name.tmp` if($txt_line > 6);
			#`cat $txt > $dir/Template/$table_name` if($txt_line <= 6);
			`head -n 37 $txt >$dir/Template/$table_name.tmp` if($txt_line > 37);
			`cat $txt > $dir/Template/$table_name` if($txt_line <= 37);
			$txt_tmp="$dir/Template/$table_name.tmp";
		}
		$txt_tmp=&T($txt_tmp);
		$content = "\t<table name=$desc type=$table_type path=$txt_tmp desc=$anno \/>\n";
=cut
		$content .="\t<file_list name=$desc1 type=$type desc=\"\" >\n";
		foreach my $path (@tables){
			my $title=&Title_write($path);
			if (!-f "$dir/HTML/$title.html"){
			  &table_html($path,"$dir/HTML/$title.html",$desc1,"","","../src/",$anno) if ($path!~/html$/ and $path!~/KEGG\.list/);
			  &table_html_cloud($path,"$dir/HTML/$title.html.cloud",$desc1,"","","project_Template_Path/../src",$anno,"logo_image_path") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
			  $path="$dir/HTML/$title.html" if ($path!~/html$/);
			}else{
			  &table_html_cloud($path,"$dir/HTML/$title.tmp.html.cloud",$desc1,"","","project_Template_Path/../src",$anno,"logo_image_path") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
			  &table_html($path,"$dir/HTML/$title.tmp.html",$desc1,"","","../src/",$anno) if ($path!~/html$/ and $path!~/KEGG\.list/);
			  $path="$dir/HTML/$title.tmp.html" if ($path!~/html$/);
			}
			my $html_title=basename($path);
			$content .="\t\t".'<file name="'.$html_title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
		}
	  $content .="\t".'</file_list>'."\n";
=cut
	}
	return $content;
  }


sub stat_table_write{
	my ($name,$type,$context) = @_;
	my ($idir,$desc,$anno);
	$desc=$$context[0];
	my $desc1=$desc;
	$desc1=~s/表\d+//;
	$desc=~s/表(\d+)\s/表$tab_nummber /  if ($desc!~/附表/);
	$tab_nummber++ if ($desc!~/附表/) ;
	$desc=&T($desc);
	$desc1=&T($desc1);
	$type = "xls";
	$type = &T($type);
	$anno=(@$context==2) ? $$context[1] :" ";
	$anno=&T($anno);
	$idir=(split /\s+/,$name)[1];
	my $file="$dir_abs/$idir";
	my $txt = (glob $file)[0];
	my $content;
	my $txt_line=`wc -l $txt|awk '{print \$1}'`; chomp $txt_line;
	if ($txt_line < 20 && $txt_line > 1) {
        $txt=&T($txt);
		$content="\t<table name=$desc type=$table_type path=$txt desc=$anno \/>\n";
    }elsif ($txt_line > 20 ){
		my $table_name=&Title_write($txt);
		my $txt_tmp=shift;
		if(!-f "$dir/Template/$table_name"){
		  `head -n 37 $txt >$dir/Template/$table_name`;
		  $txt_tmp="$dir/Template/$table_name";
		}else{
		  `head -n 37  $txt >$dir/Template/$table_name.tmp`;
		  $txt_tmp="$dir/Template/$table_name.tmp";
		}
		$txt_tmp=&T($txt_tmp);
		$content = "\t<table name=$desc type=$table_type path=$txt_tmp desc=$anno \/>\n";
		$content .="\t<file_list name=$desc1 type=$type desc=\"\" >\n";
		my $title=&Title_write($txt);
		if (!-f "$dir/HTML/$title.html"){
		  &table_html_cloud($txt,"$dir/HTML/$title.html.cloud","$desc1","","","project_Template_Path/../src",$anno,"logo_image_path") if ($txt!~/html(\.cloud)?$/ and $txt!~/KEGG\.list/);
		  &table_html($txt,"$dir/HTML/$title.html","$desc1","","","../src/",$anno) if ($txt!~/html$/ and $txt!~/KEGG\.list/);
		  $txt="$dir/HTML/$title.html" if ($txt!~/html$/);
		}else{
		  &table_html_cloud($txt,"$dir/HTML/$title.tmp.html.cloud","$desc1","","","project_Template_Path/../src",$anno,"logo_image_path") if ($txt!~/html(\.cloud)?$/ and $txt!~/KEGG\.list/);
			&table_html($txt,"$dir/HTML/$title.tmp.html","$desc1","","","../src/",$anno) if ($txt!~/html$/ and $txt!~/KEGG\.list/);
			$txt="$dir/HTML/$title.tmp.html" if ($txt!~/html$/);
		}
		my $html_title=basename($txt);
		$content .="\t\t".'<file name="'.$html_title.'" desc="" type="xls" path="'.$txt.'" action="'.'xls'.'" />'."\n";
		$content .="\t".'</file_list>'."\n";
	}
    
	
}
sub Web_Report_directory{
	my ($name,$type,$context) = @_;
	my ($idir,$desc,$anno);
	$desc=$$context[0];
	$anno=(@$context==2)?$$context[1]:" ";
	$idir=(split /\s+/,$name)[1];
	my $file="$dir_abs/$idir";
	my @tables = glob ($file);
	my $file_nu =scalar @tables;
	if ($file_nu ==0){
		print "warn $file is not exists !!";
	}
	$desc = &T($desc);
	$anno=&T($anno);
	$type=&T($type);
	my $path=&T($file);
	my $action=&T("type1");
	my $content="\t<file name=$desc type=$type path=$path action=$action />\n";
	return $content;
}
sub HTML_write{
	my ($name,$type,$context) = @_;
	my ($idir,$desc,$anno);
	$desc=$$context[0];
	$anno=(@$context==2)?$$context[1]:" ";
	$idir=(split /\s+/,$name)[1];
	my $file="$dir_abs/$idir";
	my @tables = glob ($file);
	my $file_nu =scalar @tables;
	next if $file_nu ==0;
	$desc = &T($desc);
	$anno=&T($anno);
	$type=&T($type);
	my $content="\t<file_list name=$desc type=$type desc=$anno >\n";
	foreach my $table(@tables){
		next if (!-f $table);
		my $table_name;
		if ($table=~/testForDEU.html/){
			$table=~/\/(\w+_vs_\w+)testForDEU.html/;
			$table_name=$1;
			print "$table_name\n";
			$table_name=$table_name.".DEU.html";
		}else{
			$table_name=basename($table);
		}
		$table_name=&T($table_name);
		my $path=&T($table);
		my $action=&T("type1");
		$content .= "\t\t<file name=$table_name type=$type path=$path action=$action />\n";
	}
	$content .= "\t".'</file_list>'."\n" ;
	return $content;
}

sub Title_write{
  my ($path)=@_;
  my $title=basename($path);
  if ($path=~/Cis_Anno_enrichment/){
	$title="LncRNA_Cis.".$title;
  }elsif($path=~/Trans_Anno_enrichment/){
	$title="LncRNA_Trans.".$title;
  }
  return $title;
}
#####2016-03-03
sub get_file {#
    chomp (my $id = shift);
	chomp (my $file = shift);
	if ($id=~/CellRanger分析结果网页版报告链接：/) {
		my @web_result=glob("$dir/$file");
	   $file_all_info{"CellRanger分析结果网页版报告链接："} =[@web_result];
		return $file_all_info{$id};
	}
	if ($id=~/VDJ分析结果网页版报告链接：/) {
		my @web_result=glob("$dir/$file");
		$file_all_info{"VDJ分析结果网页版报告链接："} =[@web_result];
		return $file_all_info{$id};
	}

	if ($id=~/CellRanger聚类结果：/) {
                my @web_result=glob("$dir/$file");
           $file_all_info{"CellRanger聚类结果："} =[@web_result];
                return $file_all_info{$id};
        }
	if($id=~/CDR3单链多样性分析结果/){
		my @web_result=glob("$dir/$file");
		$file_all_info{"CDR3单链多样性分析结果"} =[@web_result];
		return $file_all_info{$id};
	}
	if ($id=~/差异表达基因注释文件/) {
		 my @web_result=glob("$dir/$file");
		$file_all_info{"差异表达基因注释文件"} =[@web_result];
		return $file_all_info{$id};
	}
	if ($id=~/各样本的差异表达结果统计：/){
		my @web_result=glob("$dir/$file");
                $file_all_info{"各样本的差异表达结果统计："} =[@web_result];
                return $file_all_info{$id};

	}
}
sub web_all_write{# 
	my ($format,$text) = @_;
	my $idir=(split /\s+/,$format)[1];
	my $file="$dir_abs/$idir";
	my $context=$$text[0];
	my @context_file = (glob $file);
    if ($context=~/CellRanger分析结果网页版报告链接：|VDJ分析结果网页版报告链接：/){
        my @context_file = &get_file($context,$idir);
        my $context_file=@context_file;
        my $line= scalar @context_file;
        my ($name)= split /\s+/,$context;
        my @file_info = @{$file_all_info{$name}};
        my $content ="\t".'<file_list name="'.$context.'" type="'.'xls'.'" desc="">'."\n";
        foreach my $path (@file_info) {
            my $title=basename$path;
	    `cp $path $dir/HTML/`;
            $path="HTML/$title";
            $content .="\t\t".'<file name="'.$title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
        }
        $content .="\t".'</file_list>'."\n";
        return $content;
    }
}
sub file_all_write{
	my ($format,$text) = @_;
        my $idir=(split /\s+/,$format)[1];
        my $file="$dir_abs/$idir";
        my $context=$$text[0];
		my @context_file = &get_file($context,$idir);
		my $context_file=@context_file;
	        my $line= scalar @context_file;
		my ($name)= split /\s+/,$context;
		my @file_info = @{$file_all_info{$name}};
		my $content ="\t".'<file_list name="'.$context.'" type="'.'xls'.'" desc="">'."\n";
		foreach my $path (@file_info) {
			my $title=basename$path;
			&table_html($path,"$dir/HTML/$title.html","$name","","","../src/",$context) if ($path!~/html$/);
			$path="$dir/HTML/$title.html" if ($path!~/html$/);
	                $content .="\t\t".'<file name="'.$title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
		}
	$content .="\t".'</file_list>'."\n";
        return $content;
}

sub Head{
	my ($name,$type,$des) = @_;
	my $class;
	if ($name eq '一级标题'){
		$class='h1';
	}
	elsif ($name eq '二级标题'){
		$class = 'h2';
	}
	elsif ($name eq '三级标题') {
		$class = 'h3';
	}
	elsif ($name eq '四级标题') {
		$class = 'h4';
	}
	
	$type = &T($type);
	my $desc=$$des[0];
	$desc = &T($desc);
	my $content = "\t<$class name=$desc type=$type desc=$desc \/>\n";
	return $content;
}
sub table_html{
	my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
	my @inputs=();
	$/="\n";
	open (IN,$input)||die $!;
	$i=1;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @tmp=split/\t+/;
		if($i<=37){
			push@inputs,\@tmp;
		}else{
			last;
		}
		$i++;
	}
	$/="》";
    my $titles;
    if($input=~/Cis_Anno_enrichment/){
        $titles=basename$input;
        $titles = "lnc_cis.$titles";
    }
    elsif($input=~/Trans_Anno_enrichment/){
        $titles=basename$input;
        $titles = "lnc_trans.$titles";
    }else{
	    $titles=basename$input;
    }
    $title=~s/"//g;
	open HH,">$outHtml" or die "$!";
	print HH <<HTML;
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
	<html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
	<meta charset="UTF-8"></meta>
	<!--[if lt IE 9]>
	<script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
	<![endif]-->
	<meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
	<meta content="width=device-width, initial-scale=1" name="viewport"></meta>
	<link href="$srcPath/css/bootstrap.min.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/index.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/js/fancyBox/jquery.fancybox.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/nav.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/raxus.css" type="text/css" rel="stylesheet" />
	<script src="$srcPath/js/jquery-1.11.3.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/nav.js" type="text/javascript"></script>
	<script src="$srcPath/js/raxus-slider.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/fancyBox/jquery.fancybox.pack.js" type="text/javascript"></script>
	<script src="$srcPath/js/fancyBox/jquery.mousewheel-3.0.6.pack.js" type="text/javascript"></script>
	<script src="$srcPath/js/bootstrap.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/ready.js" type="text/javascript"></script>
	<script src="$srcPath/js/scrolltop.js" type="text/javascript"></script>
	</head>
	<body style=\"overflow:hidden\"> 
	<div class="container shadow"  id="box"><header><img src="$srcPath/images/logo.jpg" class="pull-right" /></header>
	<div role="main" ><header><h2 id="title" class="text-center">$title</h2>
	</header>
	</div>
HTML
        if($text){
		 	$text=~s/"//g;
                print HH "<div class=\"table-responsive\" id=\"textbox\"><p>$text</p></div>\n";
        }
        print HH "<div class=\"table-responsive\" id=\"box2\"><table class=\"table table-bordered table-hover table-striped\"><thead><tr class=\"bg-info\">\n";
        if($text){
                print HH <<HTML;
<script type="text/javascript">
    var textbox=\$("#textbox").height();
    var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8-textbox;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
		 var height2 = height*0.8-textbox;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML
        }
        else{
        print HH <<HTML;
<script type="text/javascript">
        var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
        var height2 = height*0.8;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML

        }

       for (my $i=0;$i<=$#{$inputs[0]};$i++){
			print HH "<th>$inputs[0][$i]</th>\n";	
		}
        print HH "</tr></thead>\n<tbody>\n";
        for (my $k=1;$k<=$#inputs ;$k++) {
                print HH "<tr>";
                
				for (my $i=0;$i<=$#{$inputs[$k]};$i++){
                        if($linkColNum){
                                my $j=$i+1;
                                if($linkColNum=~/,$j,/){
                                        print HH "<td><a href=\"$$linkHash{$titles}{$inputs[$k][$i]}\" target=\"_blank\">$inputs[$k][$i]</a></td>";
                                }
                                else{
                                        print HH "<td>$inputs[$k][$i]</td>";
                                }
                        }
                        else{
                                print HH "<td>$inputs[$k][$i]</td>";
                        }
                }
                print HH "</tr>\n";
        }
print HH <<XGL;
</tbody>
</table>
</div>
</body>
</html>
XGL
	close HH;

}
sub table_html_cloud{
    my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text,$logo_path)=@_;
    my @inputs=();
    $/="\n";
    open (IN,$input)||die $!;
    while(<IN>){
        chomp;
        next if /^\s*$/;
        my @tmp=split/\t+/;
        push@inputs,\@tmp;
    }
    $/="》";
    my $titles;
    if($input=~/Cis_Anno_enrichment/){
        $titles=basename$input;
        $titles = "lnc_cis.$titles";
    }
    elsif($input=~/Trans_Anno_enrichment/){
        $titles=basename$input;
        $titles = "lnc_trans.$titles";
    }else{
	    $titles=basename$input;
    }
    open HH,">$outHtml" or die "$!";
	$title=~s/"//g;
    print HH <<HTML;
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
    <html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
    <meta charset="UTF-8"></meta>
    <!--[if lt IE 9]>
    <script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
    <![endif]-->
    <meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
    <meta content="width=device-width, initial-scale=1" name="viewport"></meta>
    <!--added css-->
    <link rel="stylesheet" href="$srcPath/css/amend.css"/>
    <link href="$srcPath/css/bootstrap.min.css" type="text/css" rel="stylesheet" />
    <link href="$srcPath/css/index.css" type="text/css" rel="stylesheet" />
    <link href="$srcPath/js/fancyBox/jquery.fancybox.css" type="text/css" rel="stylesheet" />
    <link href="$srcPath/css/nav.css" type="text/css" rel="stylesheet" />
    <link href="$srcPath/css/raxus.css" type="text/css" rel="stylesheet" />
    <script src="$srcPath/js/jquery-1.11.3.min.js" type="text/javascript"></script>
    <script src="$srcPath/js/nav.js" type="text/javascript"></script>
    <script src="$srcPath/js/raxus-slider.min.js" type="text/javascript"></script>
    <script src="$srcPath/js/fancyBox/jquery.fancybox.pack.js" type="text/javascript"></script>
    <script src="$srcPath/js/fancyBox/jquery.mousewheel-3.0.6.pack.js" type="text/javascript"></script>
    <script src="$srcPath/js/bootstrap.min.js" type="text/javascript"></script>
    <script src="$srcPath/js/ready.js" type="text/javascript"></script>
    <script src="$srcPath/js/scrolltop.js" type="text/javascript"></script>
    </head>
    <body>
    <div class="container shadow"><header><img src="$logo_path" class="pull-right" /></header>
    <div role="main" ><header><h2 id="title" class="text-center amend-title">$title</h2>
    </header>
    </div>
HTML
        if($text){
		 	$text=~s/"//g;
                print HH "<div class=\"table-responsive\" id=\"textbox\"><p>$text</p></div>\n";
        }
        print HH "<div class=\"table-responsive\" id=\"box2\"><table class=\"table table-bordered table-hover table-striped\"><thead><tr class=\"bg-info amend\">\n";
		      if($text){
                print HH <<HTML;
<script type="text/javascript">
    var textbox=\$("#textbox").height();
    var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8-textbox;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
		 var height2 = height*0.8-textbox;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML
        }
        else{
        print HH <<HTML;
<script type="text/javascript">
        var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
        var height2 = height*0.8;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML

        }
        for (my $i=0;$i<=$#{$inputs[0]};$i++){
            print HH "<th>$inputs[0][$i]</th>\n";    
        }
        print HH "</tr></thead>\n<tbody>\n";
        for (my $k=1;$k<=$#inputs ;$k++) {
            print HH "<tr>";
            for (my $i=0;$i<=$#{$inputs[$k]};$i++){
                if($linkColNum){
                    my $j=$i+1;
                    if($linkColNum=~/,$j,/){
                        print HH "<td><a href=\"project_Template_Path/$$linkHash{$titles}{$inputs[$k][$i]}\" target=\"_blank\">$inputs[$k][$i]</a></td>";
                    }
                    else{
                        print HH "<td>$inputs[$k][$i]</td>";
                    }
                }
                else{
                    print HH "<td>$inputs[$k][$i]</td>";
                }
            }
            print HH "</tr>\n";
        }    
print HH <<XGL;
</tbody>
</table>
</div>
</body>
</html>
XGL
    close HH;

}
sub fa_html{
	my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
	open FA,"$input" ||die $!;
	my %inputs=();
	$/=">";
	while(<FA>){
		chomp;
		next if ($_=~/^\s*$/);
		my ($id,$seq)=split(/\n+/,$_,2);
		$inputs{$id}=$seq;
	}
	close FA;
	$/="》";
	open H,">$outHtml";
	print H <<HTML;
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
	<html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
	<meta charset="UTF-8"></meta>
	<!--[if lt IE 9]>
	<script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
	<![endif]-->
	<meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
	<meta content="width=device-width, initial-scale=1" name="viewport"></meta>
	<link href="$srcPath/css/bootstrap.min.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/index.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/js/fancyBox/jquery.fancybox.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/nav.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/raxus.css" type="text/css" rel="stylesheet" />
	<script src="$srcPath/js/jquery-1.11.3.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/nav.js" type="text/javascript"></script>
	<script src="$srcPath/js/raxus-slider.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/fancyBox/jquery.fancybox.pack.js" type="text/javascript"></script>
	<script src="$srcPath/js/fancyBox/jquery.mousewheel-3.0.6.pack.js" type="text/javascript"></script>
	<script src="$srcPath/js/bootstrap.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/ready.js" type="text/javascript"></script>
	<script src="$srcPath/js/scrolltop.js" type="text/javascript"></script>
	</head>
	<body>
	<div class="container shadow"><header><img src="$srcPath/images/logo.jpg" class="pull-right" /></header>
	<div role="main" ><header><h2 id="title" class="text-center">$title</h2>
	</header>
	</div>
	<div style="word-wrap:break-word;word-break:break-all">
HTML
		if($text){
			print H "<p>$text</p>\n";
		}
		for my $key(keys %inputs){
			print H ">$key<br/>\n$inputs{$key}<br/>\n";	
		}
print H <<XGL;
	</div>
	</body>
	</html>
XGL
	close H;

}
sub Text{
	my ($name,$type,$desc) = @_;
	$type = &T($type);
	$desc = &substitute($desc);
	$desc = &T($desc);
	$desc = &Trans($desc);
	my $content = "\t<p type=$type desc=$desc \/>\n";
	return $content;
}
sub Text_mul{
	my ($name,$type,$desc) = @_;
	$type = &T($type);
	my $content;
	foreach my $de (@$desc){
		$de = &substitute($de);
		$de = &T($de);
		$de = &Trans($de);
		$content .= "\t<p type=$type desc=$de \/>\n";
	}
	return $content;
}
sub Formula{
	my ($name,$type,$desc) = @_;
	$type = &T($type);
	my $content;
	foreach my $de (@$desc){
		my $text = $de;
    	$text = '<p type="type1" desc="$$'.$text.'$$" />';
    	$text ="\t".$text."\n";
		$content .= $text;
	}
	return $content;
}
sub Anno{
	my ($name,$type,$desc) =@_;
	$type = &T($type);
	$desc = &substitute($desc);
	$desc = &T($desc);
	$desc = &Trans($desc);
	my $content = "\t<p type=$type desc=$desc \/>\n";
	return $content;
}
sub href {
        my $text =shift;
        my $name = '&lt;a href=&quot;#ref'.$text.'&quot;&gt;['.$text.']&lt;/a&gt;';
        return $name;
}
sub abstract{
    my ($name,$type,$desc) = @_;
	$type = &T($type);
	my $des=$$desc[0];
    my $p_abstract = join("","&lt;p class=&quot; p-abstract&quot; &gt;",$des,"&lt;/p&gt;");
    $p_abstract = &T($p_abstract);
    my $content = "\t<report_abstract value=$p_abstract \/>\n";
	return $content;
}
sub tab2lines {
	my $text =shift;
	my $line;
	my @lines =split(/\t/,$text);
	foreach my $each (@lines) {
		$each = $each.'<br>';
		$line= $line.$each
	}
	return $line;
}
sub Reference{
	my ($name,$type,$desc) =@_;
	$type = &T($type);
	my @lines = @$desc;
	my $content = "\t<ref_list name=\"参考文献\" type=$type desc=$desc_anno >\n";
	foreach my $line (@lines) {
		my ($id,$ref_name,$ref_lin) = (split /\t/,$line)[0,1,2];
		$id = &T($id);
		$ref_name = &T($ref_name);
		$ref_lin=&T($ref_lin);
		$content = $content."\t\t<ref id=$id name=$ref_name link=$ref_lin \/>\n";
	}
	$content = $content."\t</ref_list>\n";
	return $content;
	
}
sub Attach {
	my ($name,$type,$desc) =@_;
	$type = &T($type);
	my $content = "\t<Attach name=\"附录标题\" type=$type desc=$desc_anno \/>\n";
	return $content;
}
sub T{
	my $id=shift;
	$id = '"'.$id.'"';
	return $id;
}
sub Trans{
        my $text = shift;
        my (@ref_among_word) = $text =~ /\[(\d*)\]/g;
        if (@ref_among_word != 0) {
                for (my $j=0;$j<@ref_among_word;$j++) {
                        my $ref=$ref_among_word[$j];
                        my $ref_link = &href($ref);
                        $ref= quotemeta($ref);
                        $text =~ s/\[$ref\]/$ref_link/;
                }
        }
        return $text;
}
sub substitute {
	my $text = shift;
	$text =~ s/&/&amp;/g;
	$text =~ s/'/&apos;/g;
	$text =~ s/"/&quot;/g;
	$text =~ s/</&lt;/g;
	$text =~ s/>/&gt;/g;
	
	return $text;
}
sub GetDate {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime( time() );
    return sprintf( "%4d\/%02d\/%02d", $year + 1900, $mon + 1, $day );
}
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d年%02d月%02d日", $year+1900, $mon+1, $day);
}
sub time_format {
	my $time = shift;
	my ($year,$mon,$day) = (split /\//,$time);
	return (sprintf("%4d年%02d月%02d日", $year, $mon, $day));
}
sub USAGE {
	my $usage =<< "USAGE";
	Program: $Script
	Version: $version
	Usage:
	Options:
	-indir		<dir>Web_Report
	-config			<detail.cfg>
	-template		Web_Report/template
	-xml				<file>configtest_raw.xml
	-bio			local xml
	-cloud			cloud xml
	-h			help
USAGE
	print $usage;
	exit;
}
