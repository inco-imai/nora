#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $opt_help=0;

GetOptions(
  'help'=>\$opt_help
);

my @msgs=(
  "USAGE: <this> <in.fasta>",
  "[-h show this message]"
);

if($opt_help){
  my $msg = join("\n\t",@msgs);
  printf STDERR ("%s\n",$msg);
  exit(0);
}

my $counter=0;

# fasta

my $heading = ">";

my $name = <>;
chomp $name;
$name =~ s/^$heading//;
++$counter;
my $bases = "";
my $qval = "";
while(1){
  while(my $buf=<>){
    chomp $buf;
    if($buf =~ /^$heading/){
      printf("S\t");
      printf("%s\t",$name);
      #printf("read%08d\t",$name);
      printf("%s\n",$bases);
      $name = $buf;
      $bases= "";
      $name =~ s/^$heading//;
      ++$counter;
      last;
    }
    else{
      $bases .= $buf;
    }
  }
  if(eof){
    last;
  }
}
printf("S\t");
printf("%s\t",$name);
#printf("read%08d\t",$name);
printf("%s\n",$bases);

