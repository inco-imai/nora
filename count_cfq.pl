#!/usr/bin/perl
use strict;

if(@ARGV != 1){
  die "USAGE: <this> <in.cfq>\n";
}

while(my $name=<>){
  my $bases = <>;
  my $opt = <>;
  my $mid = <>;
  chomp $mid;
  my $len = length($mid);
  my $count=0;
  for(my $i=0; $i<$len; ++$i){
    if(substr($mid,$i,1) ne 'I'){
      ++$count;
    }
  }
  print $count,"\n";
}
