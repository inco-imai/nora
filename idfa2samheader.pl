#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV != 1){
  die "USAGE: <this> <in.idfa> > <out.samheader>\n";
}

while(<>){
  chomp;
  my $name = $_;
  $name =~ s/^>//;
  my $seq = <>;
  chomp $seq;
  my $len = length($seq);
  print "\@SQ\tSN:$name\tLN:$len\n";
}

