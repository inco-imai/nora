#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $windowsize=0;
my $coverage=20;
my $opt_fastq=0;
my $opt_help=0;
my $opt_light=0;

my $estimated_genome_size = -1;

my @msgs = (
"USAGE <this> <in.fa> -g estimated_genome_size(integer>0)",
"[-w windos_size: integer>0 (0)]",
"[-c coverage: integer>0 (20)]",
"[-q: regards input as fq]",
"[-l: use less memory (but incompatible with a pipe because $0 reads <in.fa> twice)]",
"[-h: show this message]"
);

my $error_message = join "\n\t",@msgs;
#my $error_message ="USAGE <this> <in.fa> -g estimated_genome_size(integer>0)\n\t[-w windos_size(integer>0)(500)]\n\t[-c coverage(integer>0)(20)]\n\t[-q (regards input as fq)]";

GetOptions(
'g=f'=>\$estimated_genome_size,
'w=i'=>\$windowsize,
'c=i'=>\$coverage,
'q' => \$opt_fastq,
'l' => \$opt_light,
'h' => \$opt_help
);

if(@ARGV != 1){
  die "$error_message\n";
}

if($estimated_genome_size <= 0){
  die "$error_message\n";
}
if($coverage <= 0){
  die "$error_message\n";
}

my $header = ">";
if($opt_fastq){
  $header = "\@";
}

my %lengths;
my %reads;
my %others;
my %qvs;

if($opt_light){
  if($ARGV[0] eq '-'){
    printf STDERR ("-l: incompatible with a pipe because $0 reads <in.fa> twice\n");
    exit 1;
  }
  &light();
  exit 0;
}

my $buf = <>;

while(1){
  chomp $buf;
  if($buf =~ /^$header/){
    $buf =~ s/^$header//;
    my $name = $buf;
    # confirm $name was not added
    if(exists($reads{$name})){
      printf STDERR ("WARNING: the record %s is ambiguous.\n",$name);
    }
    $buf=<>;
    chomp $buf;
    my $tmp_read = $buf;
    my $n_lines = 1;
    while($buf = <>){
      chomp $buf;
      if($opt_fastq){
        if($buf =~ /^\+/){
          $reads{$name} = $tmp_read;
          last;
        }
      }
      else{
        if($buf =~ /^>/){
          $reads{$name} = $tmp_read;
          last;
        }
      }
      ++$n_lines;
      $tmp_read .= $buf;
    }
    if(eof and !$opt_fastq){
      $reads{$name} = $tmp_read;
      last;
    }
      
    if($opt_fastq){
      $others{$name} = $buf;# '+'
      my $tmp_qv = "";
      for(my $i=0; $i<$n_lines; ++$i){
        $buf = <>;
        chomp $buf;
        $tmp_qv .= $buf;
      }
      $qvs{$name} = $tmp_qv;
      if(eof){
        last;
      }
      else{
        $buf = <>;
      }
    }
  }
  else{
    printf STDERR ("strange input. $buf\n");
    exit;
  }
}

my @names = sort { length($reads{$b}) <=> length($reads{$a}) } keys %reads;

if($windowsize){
  my $maxlength=500;
  my $limit=1000000;

  for(my $i=1; $i*$windowsize <= $limit; ++$i){
    my $totalbases=0;
    my $threshold = $i*$windowsize;
    for(my $j=0; $j<@names; ++$j){
      if(length($reads{$names[$j]}) >= $threshold){
        $totalbases += length($reads{$names[$j]});
      }
      else{
        last;
      }
    }
    if($totalbases >= $coverage*$estimated_genome_size){
      $maxlength = $threshold;
    }
    else{
      last;
    }
  }

  if($maxlength+$windowsize > $limit){
    printf STDERR ("WARNING: maybe strange input.\n");
  }

  for(my $j=0; $j<@names; ++$j){
    if(length($reads{$names[$j]}) >= $maxlength){
      if(!$opt_fastq){
        printf(">%s\n",$names[$j]);
        printf("%s\n",$reads{$names[$j]});
      }
      else{
        printf("\@%s\n",$names[$j]);
        printf("%s\n",$reads{$names[$j]});
        printf("%s\n",$others{$names[$j]});
        printf("%s\n",$qvs{$names[$j]});
      }
    }
    else{
      last;
    }
  }

  printf STDERR ("%d\n",$maxlength);
}
else{
  my $total_length = 0;
  my $threshold_length = $estimated_genome_size * $coverage;
  my $prev_length = 1000000000000000;
  my $shortest_read_length=1;
  my $f_printed=0;
  for(my $j=0; $j<@names; ++$j){
    my $current_rl = length($reads{$names[$j]});
    if($current_rl > $prev_length){
      printf STDERR ("not sorted by length: %d %d\n", $prev_length, $current_rl);
      exit 1;
    }
    $prev_length = $current_rl;
    if(length($reads{$names[$j]}) >= $shortest_read_length){
      if(!$opt_fastq){
        printf(">%s\n",$names[$j]);
        printf("%s\n",$reads{$names[$j]});
      }
      else{
        printf("\@%s\n",$names[$j]);
        printf("%s\n",$reads{$names[$j]});
        printf("%s\n",$others{$names[$j]});
        printf("%s\n",$qvs{$names[$j]});
      }
      $total_length += $current_rl;
      if($total_length >= $threshold_length){
        $shortest_read_length=$current_rl;
        if(!$f_printed){
          printf STDERR ("%d\n",$current_rl);
          $f_printed = 1;
        }
        #exit 0;
      }
    }
    else{
      exit 0;
    }
  }
}

sub light(){
  open my $fh, "<$ARGV[0]" or die "cannot open $ARGV[0]: $!\n";

  my @lengths = ();
  {
    my $buf = <$fh>;
    while(1){
      chomp $buf;
      if($buf =~ /^$header/){
        $buf =~ s/^$header//;
        my $name = $buf;
        # confirm $name was not added
        if(exists($reads{$name})){
          printf STDERR ("WARNING: the record %s is ambiguous.\n",$name);
        }
        $buf=<$fh>;
        chomp $buf;
        my $tmp_read = $buf;
        my $n_lines = 1;
        while($buf = <$fh>){
          chomp $buf;
          if($opt_fastq){
            if($buf =~ /^\+/){
              #$reads{$name} = $tmp_read;
              push @lengths, length($tmp_read);
              last;
            }
          }
          else{
            if($buf =~ /^>/){
              #$reads{$name} = $tmp_read;
              push @lengths, length($tmp_read);
              last;
            }
          }
          ++$n_lines;
          $tmp_read .= $buf;
        }
        if(eof($fh) and !$opt_fastq){
          #$reads{$name} = $tmp_read;
          push @lengths, length($tmp_read);
          last;
        }

        if($opt_fastq){
          #$others{$name} = $buf;# '+'
          my $tmp_qv = "";
          for(my $i=0; $i<$n_lines; ++$i){
            $buf = <$fh>;
            chomp $buf;
            $tmp_qv .= $buf;
          }
          #$qvs{$name} = $tmp_qv;
          if(eof($fh)){
            last;
          }
          else{
            $buf = <$fh>;
          }
        }
      }
      else{
        printf STDERR ("strange input (light 1). $buf\n");
        exit 1;
      }
    }
  }
  @lengths = sort {$b <=> $a} @lengths;
  my $shortest_read_length = 1;
  my $total_length = 0;
  my $threshold_length = $estimated_genome_size * $coverage;
  my $prev_length = 1000000000000000;
  for(my $i=0; $i<@lengths; ++$i){
    my $current_rl = $lengths[$i];
    if($current_rl > $prev_length){
      printf STDERR ("not sorted by length: %d %d\n", $prev_length, $current_rl);
      exit 1;
    }
    $prev_length = $current_rl;
    $total_length += $current_rl;
    if($total_length >= $threshold_length){
      $shortest_read_length = $current_rl;
      last;
    }
  }
  if(!seek $fh,0,0){
    die "failed to seek in $0\n";
  }
  my $buf = <$fh>;
  while(1){
    my $name;
    my $bases;
    my $opts;
    my $qvs;

    chomp $buf;
    if($buf =~ /^$header/){
      $buf =~ s/^$header//;
      $name = $buf;
      $buf=<$fh>;
      chomp $buf;
      my $tmp_read = $buf;
      my $n_lines = 1;
      while($buf = <$fh>){
        chomp $buf;
        if($opt_fastq){
          if($buf =~ /^\+/){
            $bases = $tmp_read;
            last;
          }
        }
        else{
          if($buf =~ /^>/){
            $bases = $tmp_read;
            last;
          }
        }
        ++$n_lines;
        $tmp_read .= $buf;
      }
      if(eof($fh) and !$opt_fastq){
        $bases = $tmp_read;
#        &flush();
#        last;
      }

      if($opt_fastq){
        $opts = $buf; # '+'
        my $tmp_qv = "";
        for(my $i=0; $i<$n_lines; ++$i){
          $buf = <$fh>;
          chomp $buf;
          $tmp_qv .= $buf;
        }
        $qvs = $tmp_qv;
        if(eof($fh)){
#          &flush();
#          last;
        }
        else{
          $buf = <$fh>;
        }
      }
    }
    else{
      printf STDERR ("strange input (light 2). $buf\n");
      exit;
    }
    if(length($bases) >= $shortest_read_length){
      if($opt_fastq){
        printf("\@%s\n",$name);
        printf("%s\n",$bases);
        printf("%s\n",$opts);
        printf("%s\n",$qvs);
      }
      else{
        printf(">%s\n",$name);
        printf("%s\n",$bases);
      }
    }
    if(eof($fh)){
      last;
    }
  }
  close $fh;
  printf STDERR ("%d\n",$shortest_read_length);
}

