#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $prefix="";
#my $min_length=6000;
my $bin_dir="./";
#my $min_depth=5;
#my $min_depth=3;
my $min_depth=2;
#my $MIN_ALI_FRAC = 0.25;
my $MIN_ALI_FRAC = 0.0;
my $CHUNK_SIZE = 100;
my $chimeric=$CHUNK_SIZE;
#my $T_PIDEN = 0.75;
#my $T_PIDEN = 0.60;
my $T_PIDEN = 0.50;
my $genome_size = 5000000;

my $NUM_THREADS = $ENV{OMP_NUM_THREADS};

my $opt_sprai=0;

GetOptions(
  'prefix=s' =>\$prefix,
  'bin_dir=s' =>\$bin_dir,
#  'min_length=i' =>\$min_length,
  'min_depth=i' =>\$min_depth,
  't_piden=f' =>\$T_PIDEN,
  'genome_size=i' =>\$genome_size,
  'sprai' =>\$opt_sprai
);

if(@ARGV != 1 || $prefix eq ""){
  die "USAGE: <this> --prefix PREFIX <in.fa> [--genome_size int]\n";
}
my $infa=$ARGV[0];

my $ret;
my $com="";

my $inidfa="$prefix.0.idfa";
#my $out0fa="$prefix.0.done.fa";

$ret=`$bin_dir/fa2idfa_v2.pl -f --prefix $prefix.0 $infa > $inidfa`;
print "fa2idfa_v2 done\n";

$ret = `$bin_dir/get_longest_20x_fa.pl -g $genome_size -l $inidfa > $prefix.longest20x.idfa`;
print "longest20x done\n";

#core
#$ret = `$bin_dir/fa2idfa_v2.pl --prefix $prefix --valid_read_length $min_length $prefix.longest20x.idfa > $prefix.idfa`;
$ret = `$bin_dir/fa2idfa_v2.pl --prefix $prefix $prefix.longest20x.idfa > $prefix.idfa`;
print "fa2idfa_v2 done2\n";
$ret = `$bin_dir/mhseol_v3 -K 14 -H 300 -t 3 -S 5 -O -C $CHUNK_SIZE -m $MIN_ALI_FRAC -I $T_PIDEN $prefix.idfa`;
print "ol done\n";

$com="";
for(my $i=0; $i<$NUM_THREADS; ++$i){
  $com .= sprintf("$bin_dir/bfmt72s_v3 -c $chimeric -u -i -o -s $prefix.%04d.bfmt7 > $prefix.%04d.nss & ",$i,$i);
}
$ret = `$com wait`;
print "bfmt72s done\n";

$com="";
for(my $i=0; $i<$NUM_THREADS; ++$i){
  $com .= sprintf("$bin_dir/nss2v_v4 $prefix.%04d.nss > $prefix.%04d.vertical & ",$i,$i);
}
$ret = `$com wait`;
print "nss2v done\n";

$com="";
for(my $i=0; $i<$NUM_THREADS; ++$i){
  $com .= sprintf("$bin_dir/vertical2cfq_and_blacklist.pl --prefix $prefix.%04d --min_depth $min_depth --term $chimeric $prefix.%04d.vertical & ",$i,$i);
}
$ret = `$com wait`;
print "vertical2cfq_and_blacklist done\n";

$com="cat ";
for(my $i=0; $i<$NUM_THREADS; ++$i){
  $com .= sprintf("$prefix.%04d.blacklist ",$i);
}
$com .= " > $prefix.blacklist";
$ret = `$com`;
print "cat blacklists done\n";

$com="cat ";
for(my $i=0; $i<$NUM_THREADS; ++$i){
  $com .= sprintf("$prefix.%04d.cfq ",$i);
}
$com .= " > $prefix.cfq";
$ret = `$com`;
print "cat cfqs done\n";

$com="cat ";
for(my $i=0; $i<$NUM_THREADS; ++$i){
  $com .= sprintf("$prefix.%04d.bfmt7 ",$i);
}
$com .= " > $prefix.bfmt7";
$ret = `$com`;
print "cat bfmt7s done\n";

$ret = `$bin_dir/gt_v3 -B $prefix.blacklist -s $prefix -t 4 $prefix.bfmt7 > $prefix.ohf`;
print "gt done\n";
$ret = `$bin_dir/ohf2fa_v2.pl $prefix.cfq $prefix.ohf > $prefix.contigs.fa`;
$ret = `ln -s $prefix.contigs.fa $prefix.result.fa`;
print "ohf2fa done\n";
$ret = `$bin_dir/fa2idfa_v2.pl --prefix $prefix.contigs $prefix.contigs.fa > $prefix.contigs.idfa`;
print "fa2idfa done\n";
$ret = `$bin_dir/mhseol_v3 -K 14 -H 300 -t 3 -S 5 -C $CHUNK_SIZE -s $prefix.contigs.idfa > $prefix.contigs.clustering`;
print "contigs clustering done\n";

#visual
$ret = `$bin_dir/fa2gfaS.pl $prefix.idfa > $prefix.gfas`;
$ret = `cat $prefix.gfas $prefix.pruned.gfal > $prefix.pruned.gfa`;
$ret = `$bin_dir/clean_gfa.pl $prefix.pruned.gfa > $prefix.pruned.cd.gfa`;#open it using Bandage
$ret = `cat $prefix.gfas $prefix.init.gfal > $prefix.init.gfa`;
$ret = `$bin_dir/clean_gfa.pl $prefix.init.gfa > $prefix.init.cd.gfa`;
$ret = `cat $prefix.gfas $prefix.bp.gfal > $prefix.bp.gfa`;
$ret = `$bin_dir/clean_gfa.pl $prefix.bp.gfa > $prefix.bp.cd.gfa`;
$ret = `cat $prefix.gfas $prefix.orig.gfal > $prefix.orig.gfa`;
$ret = `$bin_dir/clean_gfa.pl $prefix.orig.gfa > $prefix.orig.cd.gfa`;
print "contigs gfa done\n";

#option
if($opt_sprai){
  $ret = `$bin_dir/cfq2fa.pl $prefix.cfq $prefix.blacklist > $prefix.sprai2.out.fa`;
  print "sprai2 done\n";
}

