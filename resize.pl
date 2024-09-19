#!/usr/bin/perl
use strict;

# after
#$ moddotplot static -f contigs.fa --compare-only -id 60
# , use this program.

# this program needs imagemagick
#


if(@ARGV != 2){
  die "USAGE: <this> <in.png> <out.png>\n";
}

my $in = $ARGV[0];
my $out = $ARGV[1];

my $com = "convert $in -background white -alpha off -resize 1000x1000 $out";
my $ret = `$com`;


