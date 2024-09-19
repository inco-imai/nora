#!/usr/bin/perl
use strict;

if(@ARGV != 3){
  die "USAGE: <this> <ref.fa> <query.fa> <out.blastout>\n";
}

my $ref = $ARGV[0];
my $query = $ARGV[1];
my $out = $ARGV[2];

my $command1 = "makeblastdb -in $ref -dbtype nucl -out $ref -title footitle";
`$command1`;

#my $command2 = "blastn -num_threads 6 -word_size 14 -db $ref -query $query -outfmt 5";
my $command2 = "blastn -num_threads 6 -db $ref -query $query -outfmt 5";
`$command2 > $out`;

