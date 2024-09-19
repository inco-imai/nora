#!/usr/bin/perl
use strict;
use warnings;

my $T_MERGE = 95;

if(@ARGV != 2){
  die "USAGE: <this> <contigs.fa> <get_blast_coverage.pl's output> > <result.fa>\n";
}

open my $fastafh, "<$ARGV[0]" or die "cannot open $ARGV[0] : $!\n";
open my $covfh, "<$ARGV[1]" or die "cannot open $ARGV[1] : $!\n";

my @contigs;
my $num_contigs=0;

while(my $fline = <$fastafh>){
  chomp $fline;
  $fline =~ s/^>//;
  my @tmp = split /,/,$fline;
  my $name = $tmp[0];
  #print "$fline\n$name\n";
  my $seq = <$fastafh>;
  chomp $seq;
  $contigs[$name] = $seq;
  ++$num_contigs;
}
close $fastafh;

my @scores;
for(my $i=1; $i<=$num_contigs; ++$i){
  for(my $j=1; $j<=$num_contigs; ++$j){
    $scores[$i][$j] = 0;
  }
}

while(my $cline = <$covfh>){
  chomp $cline;
  if($cline =~ /^\#/){
    next;
  }
  my @blastout = split /\t/,$cline;
  for(my $i=0; $i<=1; ++$i){
    my @tmp = split /,/,$blastout[$i];
    $blastout[$i] = $tmp[0];
  }
  if($scores[$blastout[0]][$blastout[1]] == 0){
    $scores[$blastout[0]][$blastout[1]] = $blastout[2];
  }
  elsif($scores[$blastout[0]][$blastout[1]] == $blastout[2]){
    #ok
  }
  else{
    die "unexpected scores: $scores[$blastout[0]][$blastout[1]] != $blastout[2]";
  }
}
close $covfh;

for(my $i=1; $i<=$num_contigs; ++$i){
  for(my $j=$i+1; $j<=$num_contigs; ++$j){
    if($scores[$i][$j] >= $T_MERGE && $scores[$j][$i] >= $T_MERGE){
      if(length($contigs[$i]) > length($contigs[$j])){
        print ">$i\n";
        print "$contigs[$i]\n";
      }
      elsif(length($contigs[$i]) < length($contigs[$j])){
        print ">$j\n";
        print "$contigs[$j]\n";
      }
      else{
        if($contigs[$i] le $contigs[$j]){
          print ">$i\n";
          print "$contigs[$i]\n";
        }
        else{
          print ">$j\n";
          print "$contigs[$j]\n";
        }
      }
    }
  }
}




