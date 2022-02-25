#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long::Descriptive;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $Title = "xml_add_number";
my $version="1.0.0";

my($opt,$usage)=describe_options(
'Descriptive: xml add number 
Version: 1.0.0 
Contact: linhj  <linhj\@biomarker.com.cn> 
Example: perl $Script  --xmlin configtest_raw.xml  --xmlout  configtest.xml%o',
	['xmlin|x=s',''],
	['xmlout|X=s',''],
	['help|h','print usage and exit']
);

print($usage->text),exit if(defined $opt->help);
print($usage->text),exit if(!defined $opt->xmlin && !defined $opt->xmlout);

my $xmlin=abs_path($opt->xmlin);
my $xmlout=abs_path($opt->xmlout);

#my $xml="configtest_backup.xml";
#my $out="configtest.xml";

open(IN,"$xmlin")||die "open $xmlin failed!\n";
open (OUT,">$xmlout")||die "open or creat $xmlout failed!\n";
my $i = 1 ;
my $ii;
my $j = 1 ;
my $jj;
my $k = 1 ;
my $kk;
my $l = 1 ;
my $m ;
while (<IN>) {
        chomp;
        if (/\<h1\ /){
            unless(/\<h1\ name\=\"\d+\ /){
                s/name\=\"/name\=\"$i\ /;
                s/desc\=\"/desc\=\"$i\ /;
            }
            print OUT "$_\n";
            $i++;
            $j = 1 ;
        }
        elsif(/\<h1\ name\=\"\w+/){
            s/name\=\"/name\=\"$i\ /;
            s/desc\=\"/desc\=\"$i\ /;
            print OUT "$_\n";
            $i++;
            $j = 1 ;
        }
        elsif(/\<h2\ /){
            $ii = $i - 1;
            for ($m = 1;$m<100;$m++){
                if($ii == $m ){
                    s/name\=\"/name\=\"$ii\.$j\ /;
                    s/desc\=\"/desc\=\"$ii\.$j\ /;
                    $j++;
                }
            }
            print OUT "$_\n";
            $k = 1 ;
        }
        elsif(/\<h3\ /){
            $jj = $j - 1;
            for ($m = 1;$m<100;$m++){
                if($jj == $m ){
                    s/name\=\"/name\=\"$ii\.$jj\.$k\ /;
                    s/desc\=\"/desc\=\"$ii\.$jj\.$k\ /;
                    $k++;
                }
            }
            print OUT "$_\n";
            $l = 1 ;
        }
        elsif(/\<h4\ /){
            $kk = $k - 1;
            for ($m = 1;$m<100;$m++){
                if($kk == $m ){
                    s/name\=\"/name\=\"$ii\.$jj\.$kk\.$l\ /;
                    s/desc\=\"/desc\=\"$ii\.$jj\.$kk\.$l\ /;
                    $l++;
                }
            }
            print OUT "$_\n";
        }
        else{
            #print "$_\n";
            print OUT "$_\n";
        }
    }

close IN;
close OUT;

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	use Cwd 'abs_path';
    my ($in) = @_;
    my $ret = abs_path($in);
    if ( -e $ret ){
        return $ret;
    }else{
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }
}


