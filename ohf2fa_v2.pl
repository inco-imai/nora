#!/usr/bin/perl
use strict;

if(@ARGV != 2){
  die "USAGE: <this> <in.cfq> <target.ohf>\n";
}

open (my $fhcfq, "<", $ARGV[0]) or die "cannot open $ARGV[0]\n";
open (my $fhohf, "<", $ARGV[1]) or die "cannot open $ARGV[1]\n";

#printf("%s\n",$ARGV[0]);
#printf("%s\n",$ARGV[1]);

my @acgts;
my @mids;
my $heading = "@";

my $name = <$fhcfq>;
chomp $name;
$name =~ s/^$heading//;
if($name !~ /\d+/){
  die "strange name: $name . name must be in integer.\n";
}

my $bases = "";
my $midl = "";

while(1){
  while(my $buf=<$fhcfq>){
    chomp $buf;
    #<$fhcfq>;
    #my $buf2 = <$fhcfq>;
    #chomp $buf2;
    if($buf =~ /^$heading/){
      $name =~ s/^$heading//;
      $bases =~ tr/a-z/A-Z/;
      $midl =~ tr/a-z/A-Z/;
      $acgts[$name-1] = $bases;
      $mids[$name-1] = $midl;
      if(length($acgts[$name-1]) != length($mids[$name-1])){
        printf STDERR "buggy\n";
        exit(1);
      }
      $name = $buf;
      $bases = "";
      $midl = "";
    }
    else{
      $bases .= $buf;
      <$fhcfq>;
      my $buf2 = <$fhcfq>;
      chomp $buf2;
      $midl .= $buf2;
    }
  }
  if(eof($fhcfq)){
    last;
  }
}

$name =~ s/^$heading//;
#print STDERR "#$name#\n";
$acgts[$name-1] = $bases;
$mids[$name-1] = $midl;
if(length($acgts[$name-1]) != length($mids[$name-1])){
  printf STDERR "buggy\n";
  exit(1);
}

#for(my $i=0; $i<@acgts; ++$i){
#  printf(">%d\n",$i+1);
#  printf("%s\n",$acgts[$i]);
#}

close $fhcfq;

#my $name2 = <$fhohf>;
#chomp $name2;
#$name2 =~ s/^$heading//;

my @records;
my $len = 0;

$heading = ">";

while(1){
  while(my $buf=<$fhohf>){
    chomp $buf;
    if($buf =~ /^#footer/){
      #printf("buf: %s\n",$buf);
      #printf("len %d\n",$len);
      my $name2 = $records[0];
      $name2 =~ s/^$heading//;
      my $bases2 = "";
      for(my $i=1; $i < $len; ++$i){
        #printf("records: %s\n",$records[$i]);
        my @line = split /\t/,$records[$i];
        if($line[3] eq "+"){
          my $id=$line[0];
          my $stt=$line[1];
          my $end=$line[2];
          my $newstt=0;
          my $newend=0;
          my $nbase=0;
          my $tmpi=0;
          while($nbase<=$stt-1 && $tmpi < length($mids[$id-1])){
            my $tmpmid = substr($mids[$id-1],$tmpi,1);
            if($tmpmid eq 'M'){
              ++$nbase;
            }
            elsif($tmpmid eq 'D'){
              ++$nbase;
            }
            elsif($tmpmid eq 'I'){
              #++$nbase;
            }
            ++$tmpi;
          }
          $newstt = $tmpi;
          while($nbase<=$end-1 && $tmpi < length($mids[$id-1])){
            my $tmpmid = substr($mids[$id-1],$tmpi,1);
            if($tmpmid eq 'M'){
              ++$nbase;
            }
            elsif($tmpmid eq 'D'){
              ++$nbase;
            }
            elsif($tmpmid eq 'I'){
              #++$nbase;
            }
            ++$tmpi;
          }
          $newend = $tmpi;
#print STDERR "i: $i\n";
#print STDERR "len: $len\n";
#print STDERR "newstt, newend: $newstt, $newend\n";
          $bases2 .= substr($acgts[$id-1],$newstt-1,$newend-$newstt+1);
#print STDERR $bases2,"\n";

        }
        elsif($line[3] eq "-"){
          die "strange strand '-' : $line[3]:i=$i, $buf\n";
          #my $revcomp = $acgts[$line[0]-1]; $revcomp = reverse $revcomp; $revcomp =~ tr/acgtACGT/tgcaTGCA/;
          #$bases2 .= substr($revcomp,0,$line[2]-$line[1]+1);
        }
        else{
          die "strange strand: $line[3]:i=$i, $buf\n";
        }
      }
      printf(">%s\n",$name2);
      $bases2 =~ s/-//g;
      printf("%s\n",$bases2);
      $len = 0;
      $bases2 = "";
    }
    elsif($buf =~ /^#/){
      next;
    }
    else{
      $records[$len++] = $buf;
    }
  }
  if(eof($fhohf)){
    last;
  }
}

close $fhohf;

