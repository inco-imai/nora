#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $name = "";
my $bases = "";
my $min_depth = 4;
my $min_length = 100;
my $f_wovoting = 0;

GetOptions('w'=>\$f_wovoting);

if(@ARGV != 1){
  die "USAGE: <this> <in.vertical>\n\t[-w: without voting. only chop low depth regions.]\n";
}

while(<>){
  chomp;
  if($_ eq ""){
    $bases =~ s/\-//g;
    my @tmp = split /\t/,$bases;
    my $longest = "";
    for(my $i=0; $i<@tmp; ++$i){
      if(length($longest)<length($tmp[$i])){
        $longest = $tmp[$i];
      }
    }
    $bases = $longest;
    if(length($bases)>=$min_length){
      printf(">%s\n",$name);
      printf("%s\n",$bases);
    }
    else{
      # skip
    }
    $name = "";
    $bases = "";
  }
  elsif($_ =~ /^\%/){
    $name = $_;
    $name =~ s/^\%//;
  }
  else{
    if(length($_)<$min_depth){
      $bases .= "\t";
      next;
    }
    elsif($f_wovoting){
      $bases .= substr($_,0,1);
      next;
    }
    else{
      my($nA,$nC,$nG,$nT,$nHyphen)=(0,0,0,0,0);
      $_ =~ tr/a-z/A-Z/;
      my @tmp = split //,$_;
      for(my $i=0; $i<@tmp; ++$i){
        #printf("i: %d\n",$i);
        if($tmp[$i] eq 'A'){
          ++$nA;
        }
        elsif($tmp[$i] eq 'C'){
          ++$nC;
        }
        elsif($tmp[$i] eq 'G'){
          ++$nG;
        }
        elsif($tmp[$i] eq 'T'){
          ++$nT;
        }
        elsif($tmp[$i] eq '-'){
          ++$nHyphen;
        }
        elsif($tmp[$i] eq 'N'){
          ++$nA;#XXX
        }
      }
      #printf("%d %d %d %d %d\n",$nA,$nC,$nG,$nT,$nHyphen);
      my $c_most;#character
      my $d_most;#depth
      if($nA>=$nC){
        $c_most = 'A';
        $d_most = $nA;
      }
      else{
        $c_most = 'C';
        $d_most = $nC;
      }
      if($d_most < $nG){
        $c_most = 'G';
        $d_most = $nG;
      }
      if($d_most < $nT){
        $c_most = 'T';
        $d_most = $nT;
      }
      if($d_most < $nHyphen){
        $c_most = '-';
        $d_most = $nHyphen;
      }
      $bases .= $c_most;
    }
  }
}


