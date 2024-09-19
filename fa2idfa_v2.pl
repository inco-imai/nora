#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $prefix=`date +%Y%m%d%H%M%S`;
chomp $prefix;
my $dbsuffix = ".id2n";

my $valid_read_length=1;
my $opt_fwd_only = 0;

my $min_read_length=1;
my $max_read_length=1000000000;

GetOptions(
  'prefix=s' => \$prefix,
  'valid_read_length=i'=>\$valid_read_length,
  'min_read_length=i'=>\$min_read_length,
  'max_read_length=i'=>\$max_read_length,
  'f'=>\$opt_fwd_only
);

my @msg = (
  "USAGE: <this> <in.fa> > out.idfa",
  "[--prefix prefix]",
  "[--min_read_length=i (default: $min_read_length)]",
  "[--max_read_length=i (default: $max_read_length)]",
  "[-f output forward reads only (default: false)]",
);

if(@ARGV != 1){
  my $msg = join "\n\t",@msg;
  die "$msg\n";
}

my $id_head_character=">";

my $fh;
open $fh, ">", $prefix.$dbsuffix or die "cannot open $prefix.$dbsuffix:$!\n";

my $counter=0;
my $printed_line=0;

my $line = <>;
++$counter;
while(!eof){
  chomp $line;
  my($idfaline,$nameline,$baseline);
  my($idfaline2,$nameline2);
  $idfaline = sprintf("%d\t%s",$counter,substr($line,1));
  $idfaline2 = sprintf("%d\t%s",$counter+1,substr($line,1));
  #$line = sprintf(">%d",$counter);
  $nameline = sprintf(">%d",$counter);
  if(!$opt_fwd_only){
    $nameline2 = sprintf(">%d",$counter+1);
  }
  #print $line,"\n";# name
  #++$printed_line;
  
  my $bases="";
  if(eof){
    last;
  }
  $line =<>;
  my $line_c = 1;
  chomp $line;
  while(1){# read bases
    $bases .= $line;
    if(eof){
      last;
    }
    $line = <>;
    chomp $line;
    if($line =~ /^>/){
      chomp $line;
      #$nameline = $line;#TODO
      last;
    }
    else{
      ++$line_c;
    }
  }
  $baseline = $bases;

  if(length($baseline)>=$min_read_length && length($baseline)<=$max_read_length){
    print $fh $idfaline,"\n";
    print $nameline,"\n";
    ++$printed_line;
    print $baseline,"\n";
    ++$printed_line;
    ++$counter;

    if(!$opt_fwd_only){
      print $fh $idfaline2,"\n";
      print $nameline2,"\n";
      ++$printed_line;
      $baseline =~ tr/ACGTacgt/TGCAtgca/;
      $baseline = reverse $baseline;
      print $baseline,"\n";
      ++$printed_line;
      ++$counter;
    }
  }
}

close $fh;

