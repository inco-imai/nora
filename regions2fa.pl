#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $l_margin = 50;

GetOptions(
  'l_margin=i' => \$l_margin
);

if(@ARGV != 2){
  die "USAGE: <this> <in.fwd-idfa> <in.regions> > <out.fa>\n";
}

my @stts;
my @ends;
my $prev_id = -1;
my $len = 0;

open my $regfh, "<$ARGV[1]" or die "cannot open $ARGV[1]: $!";
open my $fafh, "<$ARGV[0]" or die "cannot open $ARGV[0]: $!";


while(my $region=<$regfh>){
  chomp $region;
  my ($id,$stt,$end,$qlen) = split /\t/,$region;
  if($id != $prev_id){
    $len = 0;
    $prev_id = $id;
  }

  if($stt >= $l_margin){
    $stt += $l_margin;
  }
  if($end < $qlen-$l_margin){
    $end -= $l_margin;
  }
  $stts[$id][$len] = $stt;
  $ends[$id][$len] = $end;
  ++$len;
}

close $regfh;

while(my $id=<$fafh>){
  chomp $id;
  $id =~ s/>//;
  my $seq = <$fafh>;
  chomp $seq;
  for(my $i=0; defined($stts[$id][$i]); ++$i){
    printf(">%d/%d\n",$id,$i);
    printf("%s\n",substr($seq,$stts[$id][$i],$ends[$id][$i]-$stts[$id][$i]+1));
  }
}

close $fafh;

