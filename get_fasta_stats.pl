#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX;

my $contig=0;
my $number=0;
my $value_only=0;
my $opt_print_length=0;
my $min_len = 0;
my $genome_size = 0;
my $opt_precise=0;
my $target="";
my $opt_lens=0;

my @keys;
my @vals;

GetOptions(
  'c'=>\$contig,
  'n'=>\$number,
  'v'=>\$value_only,
  'min_len=i'=>\$min_len,
  'genome_size=i'=>\$genome_size,
  'p'=>\$opt_precise,
  't=s'=>\$target,
  'l'=>\$opt_print_length,
  'lens'=>\$opt_lens
);

if($contig){
  #$number = 1;
}

#my $error_message ="USAGE: <this> <que.fa>\n\t [-c: chop N(s) or n(s)]\n\t [-n: show the ordinal number of contigs for N50]";
my @msgs=(
  "USAGE: <this> <que.fa>",
  "[-l: print each length of contig]",
  "[-lens: input is a length value file (lens1\\nlens2\\n...)]",
  "[-min_len=i: reads shorter than this value will be ignored (default: 0)]",
  "[-genome_size=i: output NG50 (default: 0)]",
  "[-c: split input by N(s) or n(s)]",
  "[-v: show value only (without header)]",
  "[-p: show precise info.]",
  "[-t=s: fill s in a target blank]"
);
my $error_message = join("\n\t",@msgs);

if(@ARGV != 1){
  die "$error_message\n";
}

my $que=$ARGV[0];
#printf("input_file\t%s\n",$que);
push @keys,'input_name';
if(!$target){
  $target = $que;
}
push @vals,$target;

my $fh;
if($que eq "-"){
  $fh = *STDIN;
}
else{
  open $fh, "<", $que or die "cannot open $que:$!\n";
}

my %read_lengths;

if(!$opt_lens){
  my $name = <$fh>;
  chomp $name;
  my $bases = "";
  
  while(1){
    while(my $buf=<$fh>){
      chomp $buf;
      if($buf =~ /^\@/){
        die "strange line\n$buf\n";
      }
      if($buf =~ /^>/){
        $name =~ s/^>//;
        # confirm $name was not added
        if(exists($read_lengths{$name})){
          printf STDERR ("WARNING: the record %s conflicted.\n",$name);
        }

#        $bases =~ s/^N+//i;
#        $bases =~ s/N+$//i;
        if($contig){
          my @contigs = split /[N,n]+/,$bases;
          for(my $i=0; $i<@contigs; ++$i){
            my $thisname = sprintf("%s_%04d",$name,$i);
            if(length($contigs[$i])>=$min_len){
              $read_lengths{$thisname} = length($contigs[$i]);
            }
          }
        }
        else{
          if(length($bases)>=$min_len){
            $read_lengths{$name} = length($bases);
          }
          else{
            # ignore this record
          }
        }

        $name = $buf;
        $bases= "";
        last;
      }
      else{
      $bases .= $buf;
      }
    }
    if(eof){
      $name =~ s/^>//;
#      $bases =~ s/^N+//i;
#      $bases =~ s/N+$//i;
      if($contig){
        my @contigs = split /[N,n]+/,$bases;
        for(my $i=0; $i<@contigs; ++$i){
          my $thisname = sprintf("%s_%04d",$name,$i);
          if(length($contigs[$i])>=$min_len){
            $read_lengths{$thisname} = length($contigs[$i]);
          }
        }
      }
      else{
        if(length($bases)>=$min_len){
          $read_lengths{$name} = length($bases);
        }
        else{
          # ignore this record
        }
      }
      last;
    }
  }
}
else{
  for(my $i=0,my $len; $len = <$fh>; ++$i){
    chomp $len;
    $read_lengths{$i} = $len;
  }
}

close $fh;

my @names = sort { $read_lengths{$b} <=> $read_lengths{$a} } keys %read_lengths;

push @keys,'total_units';
push @vals,scalar(@names);
#printf("total_units\t%d\n",scalar(@names));

my $totalbases=0;
for(my $i=0; $i<@names; ++$i){
  $totalbases+=$read_lengths{$names[$i]};
  if($opt_print_length){
    printf("%s\t%d\n",$names[$i],$read_lengths{$names[$i]});
  }
}
#printf("total_bases\t%d\n",$totalbases);
push @keys,'total_bases';
push @vals,scalar($totalbases);

#my ($n10,$n25,$n50,$n75,$n95) =
my @nxx =
(
  ceil($totalbases*0.10),
  ceil($totalbases*0.20),
  ceil($totalbases*0.25),
  ceil($totalbases*0.30),
  ceil($totalbases*0.50),
  ceil($totalbases*0.60),
  ceil($totalbases*0.70),
  ceil($totalbases*0.80),
  ceil($totalbases*0.90),
  ceil($totalbases*0.95)
);

my @ngxx =
(
  ceil($genome_size*0.10),
  ceil($genome_size*0.20),
  ceil($genome_size*0.25),
  ceil($genome_size*0.30),
  ceil($genome_size*0.50),
  ceil($genome_size*0.60),
  ceil($genome_size*0.70),
  ceil($genome_size*0.80),
  ceil($genome_size*0.90),
  ceil($genome_size*0.95)
);

my @nxx_number=();
my @ngxx_number=();

my ($max,$min)=($read_lengths{$names[0]},$read_lengths{$names[$#names]});

for(my $i=0,$totalbases=0,my $nxxidx=0; $i<@names; ++$i){
  $totalbases+=$read_lengths{$names[$i]};
  while($totalbases >= $nxx[$nxxidx]){ # if(x){a. while(x){ a.}} = while(x){a.}
    $nxx_number[$nxxidx] = $i+1; # 1-origin
    $nxx[$nxxidx++] = $read_lengths{$names[$i]};
    last if($nxxidx >= @nxx);
  }
  last if($nxxidx >= @nxx);
}

{
  my $ngxxidx;
  for(my $i=0,$totalbases=0,$ngxxidx=0; $i<@names; ++$i){
    $totalbases+=$read_lengths{$names[$i]};
    while($totalbases >= $ngxx[$ngxxidx]){ # if(x){a. while(x){ a.}} = while(x){a.}
      $ngxx_number[$ngxxidx] = $i+1; # 1-origin
      $ngxx[$ngxxidx++] = $read_lengths{$names[$i]};
      last if($ngxxidx >= @ngxx);
    }
    last if($ngxxidx >= @ngxx);
  }
  for(;$ngxxidx < @ngxx; ++$ngxxidx){
    $ngxx[$ngxxidx] = 0;
  }
}
#printf("max\t%d\n",$max);
#printf("min\t%d\n",$min);
push @keys,'max';
push @vals,$max;
push @keys,'min';
push @vals,$min;

#my @nxxnames=("N10","N25","N50","N75","N95");
my @nxxnames=("N10","N20","N25","N30","N50","N60","N70","N80","N90","N95");
my $diff=1;
if(!$opt_precise){
  $diff = 2;
}

if(!$genome_size){
  for(my $i=0; $i<@nxx; $i+=$diff){
    push @keys,$nxxnames[$i];
    push @vals,$nxx[$i];
  }
}
else{
  #my @ngxxnames=("NG10","NG25","NG50","NG75","NG95");
  my @ngxxnames=("NG10","NG20","NG25","NG30","NG50","NG60","NG70","NG80","NG90","NG95");
  for(my $i=0; $i<@ngxx; $i+=$diff){
    push @keys,$ngxxnames[$i];
    push @vals,$ngxx[$i];
  }
}

if($number){
  for(my $i=0; $i<@nxx; ++$i){
    push @keys,"#$nxxnames[$i]";
    push @vals,$nxx_number[$i];
  }
}
#for(my $i=0; $i<@nxx; ++$i){
#  if($number){
#    printf("%s\t%d\t%d\n",$nxxnames[$i],$nxx[$i],$nxx_number[$i]);
#  }
#  else{
#    printf("%s\t%d\n",$nxxnames[$i],$nxx[$i]);
#  }
#}

#print
# TODO use join
{
  my $header="";
  for(my $i=0; $i<@keys; ++$i){
    $header .= $keys[$i];
    $header .= "\t";
  }
  chop $header;
  if(!$value_only){
    print "$header\n";
  }
}
{
  my $values="";
  for(my $i=0; $i<@vals; ++$i){
    $values .= $vals[$i];
    $values .= "\t";
  }
  chop $values;
  print "$values\n";
}
