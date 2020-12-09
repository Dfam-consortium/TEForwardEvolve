#!/bin/perl
use strict;
use Data::Dumper;

#my $evalDir = "DNATransTree-1-Tigger1-eval";
#my $evalDir = "LINETree-1-L2-eval";
#my $evalDir = "LINETree-1-L2-k80-eval";
#my $evalDir = "DNATransTree-1-Tigger1-k80-eval";
#my $evalDir = "DNATransTree-1-Tigger1-K80-auto-eval";
#my $evalDir = "DNATransTree-1-Tigger1-K80-EINS-eval";
#my $evalDir = "DNATransTree-1-Tigger1-K80-LINS-eval";
#my $evalDir = "DNATransTree-1-Tigger1-K80-GINS-eval";
my $evalDir = $ARGV[0];
#


# Files:
#
# REFERENCE MSA
# gput100-train-refmsa.fa             : Copy of the reference MSA
# gput100-train-refmsa.avgKDiv        : Linup calculated avg Kimura Divergence for reference MSA
# gput100-train-refmsa.cons.fa        : Linup calculated consensus for the ref MSA
# gput100-train-refmsa.cons.vs_seed   : NW Global alignment (Needle) vs seed sequence
# gput100-train-refmsa.cons.vs_self   : NW Global alignment (Needle) vs itself ( max score acheivable )
# <MISSING refmsa vs test> gput100-train-refmsa-vs-test.nhmmer_score
#
# PREDICTED MSA
# gput100-train-muscle.fa
# gput100-train-muscle.ama_score                 
# gput100-train-muscle.qscore_score 
# gput100-train-muscle.cons.fa
# gput100-train-muscle.cons.cm
# gput100-train-muscle.cons.cm_score           
# gput100-train-muscle.cons.nhmmer           
# gput100-train-muscle.cons.nhmmer_score
# gput100-train-muscle.hmm
# gput100-train-muscle.nhmmer
# gput100-train-muscle.nhmmer_score
# gput100-train-muscle.cons.vs_refmsacons
# gput100-train-muscle.trimmed-cons.fa
# gput100-train-muscle.trimmed-cons.cm
# gput100-train-muscle.trimmed-cons.cm_score           
# gput100-train-muscle.trimmed-cons.nhmmer           
# gput100-train-muscle.trimmed-cons.nhmmer_score
# gput100-train-muscle.trimmed.hmm
# gput100-train-muscle.trimmed.nhmmer
# gput100-train-muscle.trimmed.nhmmer_score
# gput100-train-muscle.trimmed-cons.vs_refmsacons
#
#
if ( ! -d "$evalDir/html" ) {
  mkdir("$evalDir/html");
}
opendir INDIR,"$evalDir" or die;
while ( my $entry = readdir(INDIR) ){
  if ( $entry =~ /rep-(\d+)/ ) {
    my $replicate = $1;
    print "Replicate $entry\n";
    opendir REPDIR,"$evalDir/$entry" or die;
    my %gputs = ();
    while ( my $repEntry = readdir(REPDIR) ){
      if ( $repEntry =~ /^gput(\d+)-train-.*fa/ && $repEntry !~ /[-\.]cons|seqs/ ) {
        push @{$gputs{$1}}, "$evalDir/rep-$replicate/$repEntry";
      }
    }
    closedir(REPDIR);
    foreach my $gp ( sort { $a <=> $b } keys(%gputs) ) {
      my $files  = "";
      foreach my $file ( sort { $a cmp $b } @{$gputs{$gp}} ) {
        $files .= " $file";
      }
      if ( -s "$evalDir/rep-$replicate/gput$gp-train-refiner.stk" && ! -s "$evalDir/rep-$replicate/gput$gp-train-refiner.fa" ) {
        system("/home/rhubley/projects/RepeatModeler/util/Linup -msa $evalDir/rep-$replicate/gput$gp-train-refiner.stk > $evalDir/rep-$replicate/gput$gp-train-refiner.fa" );
        $files .= " $evalDir/rep-$replicate/gput$gp-train-refiner.fa";
      }
      print "$gp\t$files\n";
      system("/home/rhubley/projects/RepeatModeler/util/viewMultipleMSA.pl $files");
      system("mv MultMSA.html $evalDir/html/rep-$replicate-gput$gp.html");
    }
  } 
}
closedir(INDIR);
 
opendir INDIR,"$evalDir" or die;
my %data = ();
while ( my $entry = readdir(INDIR) ){
  if ( $entry =~ /rep-(\d+)/ ) {
    my $replicate = $1;
    opendir REPDIR,"$evalDir/$entry" or die;
    while ( my $repEntry = readdir(REPDIR) ){
      if ( $repEntry =~ /^gput(\d+)-train-([^-\.]+)(.*)/ ) {
        my $generationsPerUnitTime = $1;
        my $method = $2;
        my $suffix = $3;
        if ( $suffix eq ".nhmmer_score" ) {
          $suffix = ".hmm.nhmmer_score";
        }
        my $refinerPadded = 0;
        if ( $method eq "refiner" && $suffix =~ /-padded(.*)/  ) {
          $suffix = $1;
          $refinerPadded = 1;
          $method = "refiner-padded";
        }
        if ( $suffix =~ /\.(ama_score)/ ) {
          my $scoresRef = readAMAScore("$evalDir/$entry/$repEntry");
          foreach my $key ( keys %{$scoresRef} ) {
            $data{$replicate}->{$generationsPerUnitTime}->{$method}->{$key} = $scoresRef->{$key};
          }
        }elsif ( $suffix =~ /\.(qscore_score)/ ) {
          my $scoresRef = readQScore("$evalDir/$entry/$repEntry");
          foreach my $key ( keys %{$scoresRef} ) {
            $data{$replicate}->{$generationsPerUnitTime}->{$method}->{$key} = $scoresRef->{$key};
          }
        }elsif ( $suffix =~ /\.nw_score/ ) {
          my $score = readNWScore("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$generationsPerUnitTime}->{$method}->{'cons.nw_score'} = $score;
        }elsif ( $suffix =~ /\.(\S+_score)/ ) {
          my $columnKey = $1;
          my $score = readScore("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$generationsPerUnitTime}->{$method}->{$columnKey} = $score;
        }elsif ( $suffix =~ /\.(avgKDiv)/ ) {
          my $avgKDiv = readAvgKDiv("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$generationsPerUnitTime}->{$method}->{'avgKDiv'} = $avgKDiv;
        }elsif ( $suffix =~ /\.vs_refmsacons/ ) {
          my ($rawScore, $refLen) = readGlobalAlign("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$generationsPerUnitTime}->{$method}->{'vs_refmsacons'} = $rawScore;
        }elsif ( $suffix =~ /\.vs_self/ ) {
          my ($rawScore, $refLen) = readGlobalAlign("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$generationsPerUnitTime}->{'refmsa'}->{'high_score'} = $rawScore;
        }

      }
    }
    closedir REPDIR;
  }
}
closedir INDIR;
print "Done reading...\n";


# Scale vs_refmsacons scores
my $high;
my $low;
foreach my $gput ( keys(%{$data{'1'}}) ) {
  foreach my $replicate ( keys(%data) ) {
    if ( ! defined $high || $high < $data{$replicate}->{$gput}->{'refmsa'}->{'high_score'} ) {
      $high = $data{$replicate}->{$gput}->{'refmsa'}->{'high_score'};
    }
    foreach my $method ( keys(%{$data{$replicate}->{$gput}}) ) {
      next if ( $method eq "refmsa" || $method eq "refiner-padded" );
      if ( ! defined $low || $low >  $data{$replicate}->{$gput}->{$method}->{'vs_refmsacons'} ) {
        $low = $data{$replicate}->{$gput}->{$method}->{'vs_refmsacons'};
      }
    }
  }
}
my $spread = abs($high - $low);
  #print "for gput $gput I found the high = $high and low = $low and spread = $spread\n";

foreach my $gput ( keys(%{$data{'1'}}) ) {
  foreach my $replicate ( keys(%data) ) {
    foreach my $method ( keys(%{$data{$replicate}->{$gput}}) ) {
       next if ( $method eq "refmsa" || $method eq "refiner-padded" );
       if ( exists $data{$replicate}->{$gput}->{$method}->{'vs_refmsacons'} ) {
         my $val = $data{$replicate}->{$gput}->{$method}->{'vs_refmsacons'};
         my $newVal = ($val - $low) / $spread;
          $data{$replicate}->{$gput}->{$method}->{'vs_refmsacons'} = $newVal;
         if ( $newVal < 0 ) { 
           print "Hmm...rep=$replicate, gput=$gput, rawscore = $val, newval = $newVal high = $high,  low=$low, spread = $spread\n";
         }
       }
    }
  }
}

# Calculate the mean and standard deviation
my %stats = ();
foreach my $replicate ( keys(%data) ) {
  foreach my $gput ( keys(%{$data{$replicate}}) ) {
    foreach my $method ( keys(%{$data{$replicate}->{$gput}}) ) {
      #next if ( $method eq "refmsa" );
      foreach my $column ( keys(%{$data{$replicate}->{$gput}->{$method}}) ) {
        $stats{$gput}->{$method}->{"sum_".$column} += $data{$replicate}->{$gput}->{$method}->{$column};
        $stats{$gput}->{$method}->{"count_".$column}++;
        $stats{$gput}->{$method}->{"mean_".$column} = $stats{$gput}->{$method}->{"sum_".$column} / $stats{$gput}->{$method}->{"count_".$column};
      }
    }
  }
}
foreach my $replicate ( keys(%data) ) {
  foreach my $gput ( keys(%{$data{$replicate}}) ) {
    foreach my $method ( keys(%{$data{$replicate}->{$gput}}) ) {
      #next if ( $method eq "refmsa" );
      foreach my $column ( keys(%{$data{$replicate}->{$gput}->{$method}}) ) {
        $stats{$gput}->{$method}->{"sum_stdev_".$column} += ($data{$replicate}->{$gput}->{$method}->{$column} - $stats{$gput}->{$method}->{"mean_".$column})**2;
        $stats{$gput}->{$method}->{"stdev_".$column} += sqrt($stats{$gput}->{$method}->{"sum_stdev_".$column}/$stats{$gput}->{$method}->{"count_".$column});
        $stats{$gput}->{$method}->{"low_stdev_".$column} = $stats{$gput}->{$method}->{"mean_".$column} - $stats{$gput}->{$method}->{"stdev_".$column};
        $stats{$gput}->{$method}->{"high_stdev_".$column} = $stats{$gput}->{$method}->{"mean_".$column} + $stats{$gput}->{$method}->{"stdev_".$column};
      }
    }
  }
}

#print "Dumper: " . Dumper(\%data) . "\n";
#print "Dumper: " . Dumper(\%stats) . "\n";
#die;


my %tables = ();
# Summary Table
foreach my $characteristic ( 'AMA_similarity_score', 'AMA_predictive_value', 
                             'QScore_Q', 'QScore_TC', 'cons.cm_score', 'cons.nhmmer_score', 'hmm.nhmmer_score', 'vs_refmsacons' ) {
  my @header = ( 'gput' );
  my @values = ();
  my $headerFlag = 1;
  foreach my $gput ( sort {$a <=> $b} keys(%stats) ) {
    my @row = ( $gput );
    foreach my $method ( 'refiner-padded', 'muscle', 'refiner', 'mafft', 'clustalw2' ) {
      next if ( ($method eq "refiner") && $characteristic !~ /(hmm\.|cons\.|vs_refmsacons)/ );
      next if ( $method eq "refiner-padded" && $characteristic eq "vs_refmsacons" );
      foreach my $stat ( 'mean', 'low_stdev', 'high_stdev' ) {
        if ( $headerFlag ) {
          push @header, "$method" . "_$stat"; 
        }
        push @row, $stats{$gput}->{$method}->{$stat . "_" . $characteristic};
      }
    }
    push @{$tables{'summary'}->{$characteristic}->{'data'}}, [@row];
    $headerFlag = 0;
  }
  push @{$tables{'summary'}->{$characteristic}->{'header'}},(@header);
}

# AvgKDiv
my @header = ( 'gput' );
my @values = ();
my $headerFlag = 1;
foreach my $gput ( sort {$a <=> $b} keys(%stats) ) {
  my @row = ( $gput );
  foreach my $stat ( 'mean', 'low_stdev', 'high_stdev' ) {
    if ( $headerFlag ) {
      push @header, "refmsa" . "_$stat"; 
    }
    push @row, $stats{$gput}->{'refmsa'}->{$stat . "_avgKDiv"};
  }
  push @{$tables{'summary'}->{'avgKDiv'}->{'data'}}, [@row];
  $headerFlag = 0;
}
push @{$tables{'summary'}->{'avgKDiv'}->{'header'}},(@header);


# Print summary table
open OUT,">$evalDir/summary.csv" or die;
foreach my $tableType ( 'summary' ) {
  foreach my $characteristic ( sort { $a cmp $b } keys %{$tables{$tableType}} ) {
    print OUT "\nTITLE: Summary of " . $characteristic . "\n";
    print OUT "" . join(",",@{$tables{$tableType}->{$characteristic}->{'header'}}) . "\n";
    foreach my $row ( @{$tables{$tableType}->{$characteristic}->{'data'}} )  {
      print OUT "" . join(",",@{$row}) . "\n";
    }
  }
}
close OUT;
  
# Replicate tables
foreach my $replicate ( sort {$a <=> $b} keys(%data) ) {
  foreach my $characteristic ( 'AMA_similarity_score', 'AMA_predictive_value',
                               'QScore_Q', 'QScore_TC', 'cons.cm_score', 'cons.nhmmer_score', 'hmm.nhmmer_score', 'vs_refmsacons' ) {
    my @header = ( 'gput' );
    my @values = ();
    my $headerFlag = 1;
    foreach my $gput ( sort {$a <=> $b} keys(%stats) ) {
      my @row = ( $gput );
      foreach my $method ( 'refiner-padded', 'muscle', 'refiner', 'mafft', 'clustalw2' ) {
        next if ( ($method eq "refiner") && $characteristic !~ /(hmm\.|cons\.|vs_refmsacons)/ );
        next if ( $method eq "refiner-padded" && $characteristic eq "vs_refmsacons" );
        if ( $headerFlag ) {
          push @header, $method
        }
        push @row, $data{$replicate}->{$gput}->{$method}->{$characteristic};
      }
      push @{$tables{$replicate}->{$characteristic}->{'data'}}, [@row];
      $headerFlag = 0;
    }
    push @{$tables{$replicate}->{$characteristic}->{'header'}},(@header);
  }
}

# avgKDiv
foreach my $replicate ( sort {$a <=> $b} keys(%data) ) {
  my @header = ( 'gput' );
  my @values = ();
  my $headerFlag = 1;
  foreach my $gput ( sort {$a <=> $b} keys(%stats) ) {
    my @row = ( $gput );
    if ( $headerFlag ) {
      push @header, 'refmsa'
    }
    push @row, $data{$replicate}->{$gput}->{'refmsa'}->{'avgKDiv'};
    push @{$tables{$replicate}->{'avgKDiv'}->{'data'}}, [@row];
    $headerFlag = 0;
  }
  push @{$tables{$replicate}->{'avgKDiv'}->{'header'}},(@header);
}

# Print replicate table
open OUT,">$evalDir/replicates.csv" or die;
foreach my $tableType ( sort {$a <=> $b} keys(%data) ) {
  my @combined_header = ();
  foreach my $characteristic ( sort {$a cmp $b} keys %{$tables{$tableType}} ) {
    foreach my $header ( @{$tables{$tableType}->{$characteristic}->{'header'}} ) {
      push @combined_header, "$characteristic:$header";
    }
  }
  print OUT "" . join(",", @combined_header) . "\n";
  my $rowIdx = 0;
  my $characteristic;
  do {
    my @combined_data = ();
    foreach $characteristic ( sort {$a cmp $b} keys %{$tables{$tableType}} ) {
      push @combined_data, @{$tables{$tableType}->{$characteristic}->{'data'}->[$rowIdx]};
    }
    print OUT "" . join(",",@combined_data) . "\n";
    $rowIdx++;
  }while ( $rowIdx <= $#{$tables{$tableType}->{'AMA_similarity_score'}->{'data'}} );
}
close OUT;


    
sub readScore {
  my $file = shift;
  
  my $score = `cat $file | tr -d '\n\r'`;
  return $score;
}
    
sub readAMAScore {
  my $file = shift;

  #Alignment 'reference.stk' has 3334 columns, 4995787 residue pairs and 3463782 insertions and 2915779 deletions
  #Alignment 'gput100-train-muscle.stk' has 3338 columns, 5037859 residue pairs and 3421710 insertions and 2873707 deletions
  #Overlap: 1325 columns, 4939231 residue pairs 3396825 insertions and 2842977 deletions
  #Normalized AMA Similarity Score:        0.984554
  #Sum of Pairs Score:                     0.988679
  #Positive Predictive Value:              0.980423
  #Total Column Score:                     0.397421
  my %scores = ();
  open IN,"<$file" or die;
  while ( <IN> ) {
    $scores{'AMA_similarity_score'} = $1 if ( /^Normalized AMA Similarity Score:\s+([\.\d]+)/ );
    $scores{'AMA_sum_of_pairs_score'} = $1 if ( /^Sum of Pairs Score:\s+([\.\d]+)/ );
    $scores{'AMA_predictive_value'} = $1 if ( /^Positive Predictive Value:\s+([\.\d]+)/ );
    $scores{'AMA_total_column_score'} = $1 if ( /^Total Column Score:\s+([\.\d]+)/ );
  }
  return \%scores;
}

sub readQScore {
  my $file = shift;

  #Test=gput100-train-muscle.fa;Ref=gput100-train-refmsa.fa;Q=0.989;TC=0.648
  my %scores = ();
  open IN,"<$file" or die;
  while ( <IN> ) {
    if ( /^Test=\S+;Ref=\S+;Q=([\.\d]+);TC=([\.\d]+)/ ) {
      $scores{'QScore_Q'} = $1;
      $scores{'QScore_TC'} = $2;
    }
  }
  return \%scores;
}

sub readNWScore {
  my $file = shift;

  # Sum Score = 
  open IN,"<$file" or die;
  my $score = "";
  while ( <IN> ) {
    if ( /^\s*Sum Score =\s+([-\d\.]+)/ ) { 
      $score = $1;
    }
  }
  return $score;
}

sub readAvgKDiv {
  my $file = shift;

  # Avg Kimura Div: 0.27
  open IN,"<$file" or die;
  my $div = "";
  while ( <IN> ) {
    if ( /^\s*Avg Kimura Div:\s+([\d\.]+)/ ) { 
      $div = $1;
    }
  }
  return $div;

}

sub readGlobalAlign {
  my $file = shift;

  # Raw Score = -188656.0
  # ..
  # reference cons length = 2337
  open IN,"<$file" or die;
  my $score = "";
  my $refLen = "";
  while ( <IN> ) {
    if ( /^\s*Raw Score =\s+([\-\d\.]+)/ ) { 
      $score = $1;
    }
    if ( /^reference cons length = (\d+)/ ) {
      $refLen = $1;
    }
  }
  return $score, $refLen;
}

