#!/usr/bin/perl
use strict;

if(@ARGV != 1){
  die "USAGE: <this> <in.cfq> > <out.fa>\n";
}

while(my $name=<>){
  chomp $name;
  my $bases = <>;
  chomp $bases;
  my $options = <>;
  chomp $options;
  my $mdis = <>;
  chomp $mdis;
  $name =~ s/^\@//;
  $bases =~ s/-//g;
  printf(">%s\n",$name);
  printf("%s\n",$bases);
}

