#!/bin/perl
use strict;
use Data::Dumper;
use JSON;


my $evalDir = $ARGV[0];

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
if ( 1 ) {
generateMSAViz( $evalDir );
}
 
opendir INDIR,"$evalDir" or die;
my %data = ();
my $xParamName;
while ( my $entry = readdir(INDIR) ){
  if ( $entry =~ /rep-(\d+)/ ) {
    my $replicate = $1;
    opendir REPDIR,"$evalDir/$entry" or die;
    while ( my $repEntry = readdir(REPDIR) ){
      if ( $repEntry =~ /^(gput|frag)(\S+)-train-([^-\.]+)(.*)/ ) {
        $xParamName = $1;
        my $xparam = $2;
        my $method = $3;
        my $suffix = $4;
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
            $data{$replicate}->{$xparam}->{$method}->{$key} = $scoresRef->{$key};
          }
        }elsif ( $suffix =~ /\.(qscore_score)/ ) {
          my $scoresRef = readQScore("$evalDir/$entry/$repEntry");
          foreach my $key ( keys %{$scoresRef} ) {
            $data{$replicate}->{$xparam}->{$method}->{$key} = $scoresRef->{$key};
          }
        }elsif ( $suffix =~ /\.nw_score/ ) {
          my $score = readNWScore("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$xparam}->{$method}->{'cons.nw_score'} = $score;
        }elsif ( $suffix =~ /\.(\S+_score)/ ) {
          my $columnKey = $1;
          my $score = readScore("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$xparam}->{$method}->{$columnKey} = $score;
        }elsif ( $suffix =~ /\.(avgKDiv)/ ) {
          my $avgKDiv = readAvgKDiv("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$xparam}->{$method}->{'avgKDiv'} = $avgKDiv;
        }elsif ( $suffix =~ /\.vs_refmsacons/ ) {
          my ($rawScore, $refLen) = readGlobalAlign("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$xparam}->{$method}->{'vs_refmsacons'} = $rawScore;
        }elsif ( $suffix =~ /\.vs_self/ ) {
          my ($rawScore, $refLen) = readGlobalAlign("$evalDir/$entry/$repEntry");
          $data{$replicate}->{$xparam}->{'refmsa'}->{'high_score'} = $rawScore;
        }

      }
    }
    closedir REPDIR;
  }
}
closedir INDIR;
print "Done reading...\n";


##
## Scale vs_refmsacons scores
##
my $high;
my $low;
# Determine Min/Max
foreach my $xval ( keys(%{$data{'1'}}) ) {
  foreach my $replicate ( keys(%data) ) {
    if ( ! defined $high || $high < $data{$replicate}->{$xval}->{'refmsa'}->{'high_score'} ) {
      $high = $data{$replicate}->{$xval}->{'refmsa'}->{'high_score'};
    }
    foreach my $method ( keys(%{$data{$replicate}->{$xval}}) ) {
      next if ( $method eq "refmsa" || $method eq "refiner-padded" );
      if ( ! defined $low || $low >  $data{$replicate}->{$xval}->{$method}->{'vs_refmsacons'} ) {
        $low = $data{$replicate}->{$xval}->{$method}->{'vs_refmsacons'};
      }
    }
  }
}
my $spread = abs($high - $low);
# Normalize values
foreach my $xval ( keys(%{$data{'1'}}) ) {
  foreach my $replicate ( keys(%data) ) {
    foreach my $method ( keys(%{$data{$replicate}->{$xval}}) ) {
       next if ( $method eq "refmsa" || $method eq "refiner-padded" );
       if ( exists $data{$replicate}->{$xval}->{$method}->{'vs_refmsacons'} ) {
         my $val = $data{$replicate}->{$xval}->{$method}->{'vs_refmsacons'};
         my $newVal = ($val - $low) / $spread;
          $data{$replicate}->{$xval}->{$method}->{'vs_refmsacons'} = $newVal;
         if ( $newVal < 0 ) { 
           print "Hmm...rep=$replicate, xval=$xval, rawscore = $val, newval = $newVal high = $high,  low=$low, spread = $spread\n";
         }
       }
    }
  }
}


##
## Calculate the mean and standard deviation
##
my %stats = ();
foreach my $replicate ( keys(%data) ) {
  foreach my $xval ( keys(%{$data{$replicate}}) ) {
    foreach my $method ( keys(%{$data{$replicate}->{$xval}}) ) {
      foreach my $column ( keys(%{$data{$replicate}->{$xval}->{$method}}) ) {
        $stats{$xval}->{$method}->{"sum_".$column} += $data{$replicate}->{$xval}->{$method}->{$column};
        $stats{$xval}->{$method}->{"count_".$column}++;
        $stats{$xval}->{$method}->{"mean_".$column} = $stats{$xval}->{$method}->{"sum_".$column} / $stats{$xval}->{$method}->{"count_".$column};
      }
    }
  }
}
foreach my $replicate ( keys(%data) ) {
  foreach my $xval ( keys(%{$data{$replicate}}) ) {
    foreach my $method ( keys(%{$data{$replicate}->{$xval}}) ) {
      foreach my $column ( keys(%{$data{$replicate}->{$xval}->{$method}}) ) {
        $stats{$xval}->{$method}->{"sum_stdev_".$column} += ($data{$replicate}->{$xval}->{$method}->{$column} - $stats{$xval}->{$method}->{"mean_".$column})**2;
        $stats{$xval}->{$method}->{"stdev_".$column} += sqrt($stats{$xval}->{$method}->{"sum_stdev_".$column}/$stats{$xval}->{$method}->{"count_".$column});
        $stats{$xval}->{$method}->{"low_stdev_".$column} = $stats{$xval}->{$method}->{"mean_".$column} - $stats{$xval}->{$method}->{"stdev_".$column};
        $stats{$xval}->{$method}->{"high_stdev_".$column} = $stats{$xval}->{$method}->{"mean_".$column} + $stats{$xval}->{$method}->{"stdev_".$column};
      }
    }
  }
}

#print "Dumper: " . Dumper(\%data) . "\n";
#print "Dumper: " . Dumper(\%stats) . "\n";
#die;


##
## Start Generating Output
##
 


my %tables = ();
# Summary Table
foreach my $characteristic ( 'AMA_similarity_score', 'AMA_predictive_value', 
                             'QScore_Q', 'QScore_TC', 'cons.cm_score', 'cons.nhmmer_score', 'hmm.nhmmer_score', 'vs_refmsacons' ) {
  my @header = ( $xParamName );
  my @values = ();
  my $headerFlag = 1;
  foreach my $xval ( sort {$a <=> $b} keys(%stats) ) {
    my @row = ( $xval );
    foreach my $method ( 'refiner-padded', 'muscle', 'refiner', 'mafft', 'clustalw2', 'dialign', 'kalign', 'fsa' ) {
      next if ( ($method eq "refiner") && $characteristic !~ /(hmm\.|cons\.|vs_refmsacons)/ );
      next if ( $method eq "refiner-padded" && $characteristic eq "vs_refmsacons" );
      foreach my $stat ( 'mean', 'low_stdev', 'high_stdev' ) {
        if ( $headerFlag ) {
          push @header, "$method" . "_$stat"; 
        }
        push @row, $stats{$xval}->{$method}->{$stat . "_" . $characteristic};
      }
    }
    push @{$tables{'summary'}->{$characteristic}->{'data'}}, [@row];
    $headerFlag = 0;
  }
  push @{$tables{'summary'}->{$characteristic}->{'header'}},(@header);
}

# AvgKDiv
my @header = ( $xParamName );
my @values = ();
my $headerFlag = 1;
foreach my $xval ( sort {$a <=> $b} keys(%stats) ) {
  my @row = ( $xval );
  foreach my $stat ( 'mean', 'low_stdev', 'high_stdev' ) {
    if ( $headerFlag ) {
      push @header, "refmsa" . "_$stat"; 
    }
    push @row, $stats{$xval}->{'refmsa'}->{$stat . "_avgKDiv"};
  }
  push @{$tables{'summary'}->{'avgKDiv'}->{'data'}}, [@row];
  $headerFlag = 0;
}
push @{$tables{'summary'}->{'avgKDiv'}->{'header'}},(@header);


# Print summary tables
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

my $xAxis = "Average Kimura Divergence";
if ( $xParamName eq "frag" ) {
  $xAxis = "Sequence Fragment Size Distribution (bp_mean - bp_stdev)";
}
my %methodColors = ( 
      'refiner-padded_mean' => "#3366CC",
      'refiner_mean' => "#0099C6",
      'muscle_mean' => "#DC3912",
      'mafft_mean' => "#FF9900",
      'clustalw2_mean' => "#109618",
      'dialign_mean' => "#990099",
      'kalign_mean' => "#AAAA11",
      'fsa_mean' => "#DD4477",
      'refmsa_mean' => "#0A69A2" ,
      'refmsa_low_stdev' => "#81B7D8",
      'refmsa_high_stdev' => "#81B7D8" );

my %summaryGraphMetaData = ( 'avgKDiv' => 
                                 { 'title' => 'Reference MSA Average Kimura Divergence',
                                   'xAxis' => 'Generations Per Unit Time',
                                   'yAxis' => 'Average Kimura Divergence' },
                             'AMA_similarity_score' => 
                                 { 'title' => 'AMA Similarity Score',
                                   'xAxis' => $xAxis,
                                   'yAxis' => 'Similarity Score' },
                             'AMA_predictive_value' => 
                                 { 'title' => 'AMA Predictive Value [modeler score]',
                                   'xAxis' => $xAxis,
                                   'yAxis' => 'Predictive Value' },
                             'QScore_Q' =>
                                 { 'title' => 'QScore Developer Score [sum of pairs]',
                                   'xAxis' => $xAxis,
                                   'yAxis' => 'Developer Score' },
                             'QScore_TC' =>
                                 { 'title' => 'QScore Total Column Score',
                                   'xAxis' => $xAxis,
                                   'yAxis' => 'Total Column Score' },
                             'cons.cm_score' =>
                                 { 'title' => 'Crossmatch - Predicted Consensus vs Test Sequences',
                                   'xAxis' => $xAxis,
                                   'yAxis' => 'Total Complexity Adjusted Score' },
                             'cons.nhmmer_score' =>
                                 { 'title' => 'Nhmmer - Predicted Consensus vs Test Sequences',
                                   'xAxis' => $xAxis,
                                   'yAxis' => 'Total Bit Score' },
                             'hmm.nhmmer_score' => 
                                 { 'title' => 'Nhmmer - Predicted HMM vs Test Sequences',
                                   'xAxis' => $xAxis,
                                   'yAxis' => 'Total Bit Score' },
                             'vs_refmsacons' =>
                                 { 'title' => 'Predicted vs Reference MSA Consensus',
                                   'xAxis' => $xAxis,
                                   'yAxis' => 'Ideal NW Score Fraction' },
                           );

# Print summary graphs
open OUT,">$evalDir/graphs.html" or die;

  print OUT "<html>\n";
  if ( $xParamName eq "frag" ) {
    print OUT "<h2>Fragmentation Evaluation</h2>\n";
  }else {
    print OUT "<h2>Divergence Evaluation</h2>\n";
  }
  print OUT "<h3>$evalDir</h3>\n";
  print OUT "<head>\n";
  print OUT "  <script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n";
  print OUT "  <script type=\"text/javascript\">\n";
  print OUT "     google.charts.load('current', {'packages':['corechart']});\n";
  print OUT "     google.charts.setOnLoadCallback(drawChart);\n";
  print OUT "     function drawChart() {\n";
  my @divs = ();
  foreach my $characteristic ( sort { $a cmp $b } keys %{$tables{'summary'}} ) {
    my $metaData = $summaryGraphMetaData{$characteristic};
    my @data = ();
    my @colors = ();
    my $dataTable = [];
    my $useStats = undef;
    if ( $characteristic =~ /KDiv/ ) {
      next if ( $xParamName eq "frag" );
      $useStats = 1;
    }
    for ( my $i = 0; $i <= $#{$tables{'summary'}->{$characteristic}->{'header'}}; $i++ ) {
      my $heading = $tables{'summary'}->{$characteristic}->{'header'}->[$i];
      next if ( !defined $useStats && $heading =~ /_stdev/ );
      push @data, $heading;
      if ( exists $methodColors{$heading} ) {
        push @colors, $methodColors{$heading};
      }
    }
    push @{$dataTable}, [@data];
    foreach my $row ( @{$tables{'summary'}->{$characteristic}->{'data'}} )  {
      @data = ();
      for ( my $i = 0; $i <= $#{$tables{'summary'}->{$characteristic}->{'header'}}; $i++ ) {
        my $heading = $tables{'summary'}->{$characteristic}->{'header'}->[$i];
        next if ( !defined $useStats && $heading =~ /_stdev/ );

        # Convert GPUT to average divergence
        my $val = $row->[$i];
        if ( $heading eq "gput" && $characteristic !~ /KDiv/ ) {
          $val =  $stats{$val}->{'refmsa'}->{"mean_avgKDiv"};
          push @data,  $val;
        }elsif ( $heading eq "frag" ) {
          push @data,  $val . "";
        }elsif ( defined $val ) {
           push @data,  $val + 0;
        }else {
           push @data, 0;
        }
      }
      push @{$dataTable},[@data];
    }

    my $graphDir = 1;
    $graphDir = -1 if ( $xParamName eq "frag" );
    my ($hashID, $jStr, $hStr) = generateChart($metaData->{'title'},$metaData->{'xAxis'}, $metaData->{'yAxis'}, $dataTable, $graphDir, \@colors);
    print OUT "        $jStr\n";
    print OUT "     var chart_$hashID = new google.visualization.LineChart(document.getElementById('div_$hashID'));\n";
    print OUT "     chart_$hashID.draw(data_$hashID,options_$hashID);\n\n";
    push @divs, $hStr;
  }
  print OUT "     }\n";
  print OUT "  </script>\n";
  print OUT "</head>\n";
  print OUT "<body>\n";
  foreach my $div ( @divs ) {
    print OUT "   $div\n";
  }
  print OUT "</body>\n";
  print OUT "</html>\n"; 
close OUT;

exit;
  
# Replicate tables
foreach my $replicate ( sort {$a <=> $b} keys(%data) ) {
  foreach my $characteristic ( 'AMA_similarity_score', 'AMA_predictive_value',
                               'QScore_Q', 'QScore_TC', 'cons.cm_score', 'cons.nhmmer_score', 'hmm.nhmmer_score', 'vs_refmsacons' ) {
    my @header = ( 'gput' );
    my @values = ();
    my $headerFlag = 1;
    foreach my $gput ( sort {$a <=> $b} keys(%stats) ) {
      my @row = ( $gput );
      foreach my $method ( 'refiner-padded', 'muscle', 'refiner', 'mafft', 'clustalw2', 'dialign', 'kalign', 'fsa' ) {
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

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################


sub generateMSAViz {
  my $evalDir = shift;

  print "Generating HTML Visualization...\n";
  if ( ! -d "$evalDir/html" ) {
    mkdir("$evalDir/html");
  }
  opendir INDIR,"$evalDir" or die;
  while ( my $entry = readdir(INDIR) ){
    if ( $entry =~ /rep-(\d+)/ ) {
      my $replicate = $1;
      print "Replicate $entry\n";
      opendir REPDIR,"$evalDir/$entry" or die;
      my %xparam = ();
      my $prefix;
      while ( my $repEntry = readdir(REPDIR) ){
        if ( $repEntry =~ /^(frag|gput)(\S+)-train-.*fa/ && $repEntry !~ /[-\.]cons|seqs/ ) {
          $prefix = $1;
          push @{$xparam{$2}}, "$evalDir/rep-$replicate/$repEntry";
        }
      }
      closedir(REPDIR);
      foreach my $gp ( sort { $a <=> $b } keys(%xparam) ) {
        my $files  = "";
        foreach my $file ( sort { $a cmp $b } @{$xparam{$gp}} ) {
          $files .= " $file";
        }
        if ( -s "$evalDir/rep-$replicate/$prefix$gp-train-refiner.stk" && ! -s "$evalDir/rep-$replicate/$prefix$gp-train-refiner.fa" ) {
          system("/home/rhubley/projects/RepeatModeler/util/Linup -msa $evalDir/rep-$replicate/$prefix$gp-train-refiner.stk > $evalDir/rep-$replicate/$prefix$gp-train-refiner.fa" );
          $files .= " $evalDir/rep-$replicate/$prefix$gp-train-refiner.fa";
        }
        print "$gp\t$files\n";
        system("/home/rhubley/projects/RepeatModeler/util/viewMultipleMSA.pl $files");
        system("mv MultMSA.html $evalDir/html/rep-$replicate-$prefix$gp.html");
      }
    } 
  }
  closedir(INDIR);
}



sub generateChart {
  my $title = shift;
  my $xTitle = shift;
  my $yTitle = shift;
  my $dataTable = shift;
  my $xDirection = shift;
  my $colors = shift;

  my $drawChartJSON = "";
  my $bodyStr = "";

  my $hashID = "";
  my @chars = ('0'..'9', 'A'..'F');
  my $len = 5;
  while($len--){ $hashID .= $chars[rand @chars] };

  $drawChartJSON .= "var data_$hashID = google.visualization.arrayToDataTable(";
  #my $data = Dumper($dataTable);
  #$data =~ s/^\$VAR1 = //;
  #$data =~ s/(.*);$/$1/g;
  my $data = encode_json($dataTable);
  $drawChartJSON .= $data;
  $drawChartJSON .= ");\n";
  $drawChartJSON .= "var options_$hashID = {\n";
  $drawChartJSON .= "  title: \'$title\',\n";
  $drawChartJSON .= "  width: 1200,\n";
  $drawChartJSON .= "  height: 500,\n";
  $drawChartJSON .= "  titleTextStyle: { color: '#808080', fontSize: 20, bold: false },\n";
  $drawChartJSON .= "  legend: { textStyle: { color: '#505050', fontSize: 12, italic: false }},\n";
  if ( $xDirection == -1 ) {
    $drawChartJSON .= "  hAxis: { title: \'$xTitle\', textStyle: { fontSize: 12 }, titleTextStyle: { color: '#505050', fontSize: 12, italic: false }, direction: -1 },\n";
  }else {
    $drawChartJSON .= "  hAxis: { title: \'$xTitle\', textStyle: { fontSize: 12 }, titleTextStyle: { color: '#505050', fontSize: 12, italic: false }},\n";
  }
  $drawChartJSON .= "  vAxis: { title: \'$yTitle\', textStyle: { fontSize: 12 }, titleTextStyle: { color: '#505050', fontSize: 12, italic: false } },\n";
  # 
  # colors: ['#AB0D06', '#007329'],... add more than enough
  $drawChartJSON .= "    colors: " . encode_json($colors) . "\n";
  $drawChartJSON .= "};\n";

  $bodyStr .= "<div id=\"div_$hashID\"></div>";

  return ($hashID, $drawChartJSON, $bodyStr);
}


    
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

