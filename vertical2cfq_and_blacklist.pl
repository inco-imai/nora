#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $name = "";
my $bases = "";
my $depths = "";
my $qvs="";
my $min_depth = 4;
my $min_length = 100;
my $f_wovoting = 0;
my $term = 300;
my $f_old_algo=0;

my $prefix = "Sprai2.out";

GetOptions(
  'w'=>\$f_wovoting,
  'prefix=s' => \$prefix,
  'min_depth=i' => \$min_depth,
  'term=i' => \$term
  );

if(@ARGV != 1){
  die "USAGE: <this> --prefix <prefix> <in.vertical>\n\t[-w: without voting. only chop low depth regions.]\n\t[-min_depth integer]\n\t[-term integer]\n";
}

my $corrected_reads = "$prefix.cfq";
my $blacklist = "$prefix.blacklist";

open my $CORRECTED, ">$corrected_reads" or die "cannot open $corrected_reads:$!\n";
open my $BLACKLIST, ">$blacklist" or die "cannot open $blacklist:$!\n";

while(<>){
  chomp;
  if($_ eq ""){
    my $stt = $term;
    my $end = length($depths)-1-$term;
    my $tmpstr = substr($depths,$stt,$end-$stt+1);
    #print STDERR $tmpstr,"\n";
    if($tmpstr =~ /,/){
      printf $BLACKLIST ("%s\n",$name);
    }
    #for(my $i=$stt; $i<=$end; ++$i){
    #  if(substr($depths,$i,1) eq ','){
    #    printf $BLACKLIST ("%s\n",$name);
    #    last;
    #  }
    #}


    if($f_old_algo){
      my @raw = split /,/,$depths;
      #my @mid = split /,/,$qvs;

      $depths =~ s/\-//g;
      my @tmp = split /,/,$depths;
      my $longest = "";
      my $second_longest = "";
      my $li=0;
      for(my $i=0; $i<@tmp; ++$i){
        if(length($longest)<length($tmp[$i])){
          $longest = $tmp[$i];
          $li = $i;
        }
      }
      for(my $i=0; $i<@tmp; ++$i){
        if($i == $li){
          next;
        }
        if(length($second_longest)<length($tmp[$i])){
          $second_longest = $tmp[$i];
        }
      }
      if(length($second_longest)>=$min_length){
        # number of valid reads >= 2
        # write the id into a blacklist
        printf $BLACKLIST ("%s\n",$name);
      }
      elsif(length($longest)>=$min_length){
        # number of valid reads == 1
        # NOT write into a blacklist
      }
      else{
        # number of valid reads == 0
        # write the id into a blacklist
        printf $BLACKLIST ("%s\n",$name);
      }
    }

    # output a cfq
    printf $CORRECTED ("@%s\n",$name);
    printf $CORRECTED ("%s\n",$bases);
    printf $CORRECTED ("+\t%s\n",$depths);
    printf $CORRECTED ("%s\n",$qvs);

    $name = "";
    $bases = "";
    $depths = "";
    $qvs = "";
  }
  elsif($_ =~ /^\%/){
    $name = $_;
    $name =~ s/^\%//;
  }
  else{
    my $tmptmp = $_;
    $tmptmp =~ s/\s//g;
    if(length($tmptmp)<$min_depth){
      $depths .= ',';# ',' means invalid in this context
      my $tmpb = substr($_,0,1);
      if($tmpb ne '-'){
        $bases .= $tmpb; #don't build a consensus
        $qvs .= 'M';
      }
      next;
    }
    else{
      my($nA,$nC,$nG,$nT,$nHyphen)=(0,0,0,0,0);
      $_ =~ tr/a-z/A-Z/;
      my @tmp = split //,$_;
      for(my $i=0; $i<@tmp; ++$i){
        #printf("i: %d\n",$i);
        if($tmp[$i] eq 'A'){
          ++$nA;
        }
        elsif($tmp[$i] eq 'C'){
          ++$nC;
        }
        elsif($tmp[$i] eq 'G'){
          ++$nG;
        }
        elsif($tmp[$i] eq 'T'){
          ++$nT;
        }
        elsif($tmp[$i] eq '-'){
          ++$nHyphen;
        }
        elsif($tmp[$i] eq 'N'){
          ++$nA;#XXX
        }
      }
      #printf("%d %d %d %d %d\n",$nA,$nC,$nG,$nT,$nHyphen);
      my $c_most;#character
      my $d_most;#depth
      if($nA>=$nC){
        $c_most = 'A';
        $d_most = $nA;
      }
      else{
        $c_most = 'C';
        $d_most = $nC;
      }
      if($d_most < $nG){
        $c_most = 'G';
        $d_most = $nG;
      }
      if($d_most < $nT){
        $c_most = 'T';
        $d_most = $nT;
      }
      if($d_most < $nHyphen){
        $c_most = '-';
        $d_most = $nHyphen;
      }
      if($f_wovoting){
        # XXX
        $bases .= $tmp[0];
        $depths .= $tmp[0];
      }
      else{
        if($tmp[0] ne '-' || $c_most ne '-'){
          $bases .= $c_most;
          $depths .= $c_most;
        }
        # don't append '-' if $tmp[0] and $c_most eq '-'
      }
      if($tmp[0] eq '-'){
        if($c_most ne '-'){
          $qvs .= 'I';
        }
      }
      elsif($c_most eq '-'){
        $qvs .= 'D';
      }
      elsif($tmp[0] eq $c_most){
        $qvs .= 'M';
      }
      else{
        $qvs .= 'm';
      }
    }
  }
}

close $CORRECTED;
close $BLACKLIST;

