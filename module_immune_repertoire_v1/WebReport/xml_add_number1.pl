#!/usr/bin/perl -w
use strict;
#use newPerlBase;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $Title = "xml_add_number";
my $version="1.0.0";

my ( $xmlin, $xmlout);
GetOptions(
    "xmlin:s"            => \$xmlin,
    "xmlout:s"            => \$xmlout,
    "help|h"          => \&USAGE,
) or &USAGE;
&USAGE unless ( $xmlin and $xmlout);
$xmlin = &ABSOLUTE_DIR($xmlin);

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

###########################################
###########################################
sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: 1.0.0
   Contact: linhj  <linhj\@biomarker.com.cn>
     Usage:
            --xmlin              input xml file              [forced ]
            --xmlout             output xml file             [forced ]
            --h                  help documents

   Example:
            perl $Script  --xmlin configtest_raw.xml  --xmlout  configtest.xml

----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}

