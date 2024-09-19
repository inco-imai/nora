#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX 'ceil';

my $line_length = 100;

GetOptions(
  'l=i'=>\$line_length
);

if(@ARGV != 1){
  die "USAGE: <this> <in.bfmt7>\n\t[-l=d: line length]\n";
}

open my $fh, "<$ARGV[0]" or die "cannot open $ARGV[0] : $!\n";

while(my $line = <$fh>){
  chomp $line;
  if($line =~ /^#/){
    next;
  }
  my @el = split /\t/,$line;
  my $qseq = $el[9];
  my $sseq = $el[10];
  my $qid = $el[0];
  my $sid = $el[3];
  print "qid>$qid\n";
  print "sid>$sid\n";
  #print "$qseq\n";
  #print "$sseq\n";
  #exit;
  my $num_l = ceil ((length($qseq)+0.0)/($line_length+0.0));
  my $n_padding = $line_length-length($qseq)%$line_length;
  for(my $i=0; $i<$n_padding; ++$i){
    $qseq .= " ";
    $sseq .= " ";
  }
  for(my $i=0; $i<$num_l; ++$i){
    my $p_q = substr($qseq,$line_length*$i,$line_length);
    my $p_s = substr($sseq,$line_length*$i,$line_length);
    my $p_a = "";
    my $piden;
    my $n_match=0;
    my $n_bases_and_indel=0;
    for(my $j=0; $j<$line_length; ++$j){
      my $e_q = substr($p_q,$j,1);
      my $e_s = substr($p_s,$j,1);
      if($e_q ne " "){
        ++$n_bases_and_indel;
        if($e_q eq $e_s){
          $p_a .= "|";
          ++$n_match;
        }
        else{
          $p_a .= " ";
        }
      }
      else{
        $p_a .= " ";
      }
    }
    printf("%d\n",$line_length*$i+1);
    print "$p_q\n";
    printf("%s\t%0.2f\n",$p_a,($n_match+0.0)/($n_bases_and_indel+0.0));
    print "$p_s\n";
    print "\n";
  }
}

