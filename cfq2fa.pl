#!/usr/bin/perl
use strict;

if(@ARGV != 2){
  die "USAGE: <this> <in.cfq> <in.blacklist> > <out.fa>\n";
}

open my $incfq, "<$ARGV[0]" or die "cannot open $ARGV[0]: $!\n";
open my $inbl, "<$ARGV[1]" or die "cannot open $ARGV[1]: $!\n";

my @blist;
my $listlen=0;
while(my $el=<$inbl>){
  chomp $el;
  $blist[$listlen] = $el;
  ++$listlen;
}
close $inbl;

while(my $name=<$incfq>){
  chomp $name;
  my $bases = <$incfq>;
  chomp $bases;
  my $options = <$incfq>;
  chomp $options;
  my $mdis = <$incfq>;
  chomp $mdis;
  $name =~ s/^\@//;
  $bases =~ s/-//g;
  if(grep {$_ eq $name} @blist){
  }
  else{
    printf(">%s\n",$name);
    printf("%s\n",$bases);
  }
}
close $incfq;

