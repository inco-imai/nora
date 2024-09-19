#!/usr/bin/perl
use strict;

my @es;
my %el;
my @others;
my $escounter = 0;
my $elcounter = 0;

open (my $fh, "<", $ARGV[0]) or die "cannot open $ARGV[0]";

while(<$fh>){
  chomp;
  my @input = split /\t/,$_;
  if ($input[0] eq 'S'){
    $es[$escounter] = $input[1];
    ++$escounter;
  }
  elsif ($input[0] eq 'L' && $input[5] ne "25M"){
    $el{$input[1]} = 1;
    $el{$input[3]} = 1;
  }
  else{
    # do nothing
  }
}

seek($fh,0,0);

while(<$fh>){
  chomp;
  my @input = split /\t/,$_;
  if ($input[0] eq 'S'){
    if ($el{$input[1]} == 1){
      printf("%s\n",$_);
    }
  }
  else{
    printf("%s\n",$_);
  }
}

