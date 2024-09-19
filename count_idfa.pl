#!/usr/bin/perl
use strict;

if(@ARGV != 1){
  die "USAGE: <this> <in.idfa>\n";
}

while(my $name=<>){
  my $bases = <>;
  chomp $bases;
  my $len = length($bases);
  print $len,"\n";
}
