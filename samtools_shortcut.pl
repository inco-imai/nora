#!/usr/bin/perl

if(@ARGV != 1){
  die "USAGE: <this> <in.sam>\n";
}

my $insam = $ARGV[0];
$insam =~ /^(\w+.*)\.sam$/;
my $prefix = $1;
#print "$prefix\n";
#exit;

my $com1 = "samtools view -bS $insam > $prefix.bam";
my $out = `$com1`;
my $com2 = "samtools sort $prefix.bam > $prefix.sorted.bam";
$out = `$com2`;
my $com3 = "samtools index $prefix.sorted.bam";
$out = `$com3`;
#my $com4 = "samtools tview $prefix.sorted.bam";
#`$com4`;

