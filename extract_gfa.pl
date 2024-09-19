#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV != 3){
  die "USAGE: <this> <in.final.fa> <in.ohf> <in.gfa> > <out.gfa>\n";
}

open my $fafh, "<$ARGV[0]" or die "cannot open $ARGV[0] : $!\n";
open my $ohffh, "<$ARGV[1]" or die "cannot open $ARGV[1] : $!\n";
open my $gfafh, "<$ARGV[2]" or die "cannot open $ARGV[2] : $!\n";

my @names;
my $num_fasta_records=0;
while(my $fl = <$fafh>){
  chomp $fl;
  $fl =~ s/^>//;
  $names[$num_fasta_records] = $fl;
  ++$num_fasta_records;
  $fl = <$fafh>; #discard
}
close $fafh;

my @records;
my $len=0;
my %exist;
while(1){
  while(my $buf=<$ohffh>){
    chomp $buf;
    if($buf =~ /^#footer/){
      my $name = $records[0];
      $name =~ s/^>//;
      my @tmp = split /,/,$name;
      $name = $tmp[0];
      #print STDERR "$name\n";
      my $flag=0;
      for(my $i=0; $i<$num_fasta_records; ++$i){
        if($names[$i] eq $name){
          $flag=1;
          last;
        }
      }
      if($flag){
        for(my $i=0; $i<$len; ++$i){
          my @tmp = split /\t/,$records[$i];
          #print STDERR "$tmp[0]\n";
          $exist{$tmp[0]} = 1;
        }
      }
      $len = 0;
    }
    elsif($buf =~ /^#/){
      next;
    }
    else{
      $records[$len++] = $buf;
    }
  }
  if(eof($ohffh)){
    last;
  }
}
close $ohffh;

while(my $gl=<$gfafh>){
  chomp $gl;
  if($gl =~ /^S/){
    print "$gl\n";
  }
  elsif($gl =~ /^L/){
    my @tmp = split /\t/,$gl;
    my $from = $tmp[1];
    my $to = $tmp[3];
    if(defined($exist{$from}) || defined($exist{$to})){
      print "$gl\n";
    }
  }
  else{
    print "$gl\n";
  }
}
close $gfafh;


