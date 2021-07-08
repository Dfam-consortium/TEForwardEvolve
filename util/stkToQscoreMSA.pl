#!/usr/bin/perl
use strict;
use Data::Dumper;

#./stkToQscoreMSA.pl reference-sequences.fa refiner-msa.stk


#
# Read the reference sequences:
#   1. Store the sequences in a hash with the identifier as the key
#   2. Initialize the missingIDs hash for later identification of missing sequences
#
my $id;
my $seq;
my %seqs;
open FA,"<$ARGV[0]" or die;
my %missingSeqIDs = ();
while ( <FA> ) {
  if ( /^>(\S+)/ ) {
    my $tID = $1;
    if ( $seq ) {
      $seqs{$id} = $seq;
      $missingSeqIDs{$id} = 1;
    }
    $id = $tID;
    $seq = "";
    next;
  }
  s/[\n\r\s]+//g;
  $seq .= $_;
}
if ( $seq ) {
  $seqs{$id} = $seq;
  $missingSeqIDs{$id} = 1;
}
close FA;

#
# Parse the multiple alignment and look for duplicate fragments or reverse-complement matches
#
open STK,"<$ARGV[1]" or die;
my %msa = ();
while ( <STK> ) {
  next if ( /^ STOCKHOLM/ );
  next if ( /^\#/ );
  next if ( /^\/\// );
  if ( /^(\S+)\s+([ACGTNRYMKSWHBVD\.]+)/ ) {
    my $originalID = $1;
    my $id = $1;
    my $seq = uc($2);
    my $start;
    my $end;
    my $aStart;
    my $aEnd;
    if ( $id =~ /^(\S+)\:(\d+)\-(\d+)/ ) {
      # Strip off additional sequence ranges to get the original identifier
      $id = $1;
      $start = $2;
      $end = $3;
    }
    # Sequence is aligned in the "wrong" orientation.
    next if ( $end < $start );

    # capture alignment start/end
    if ( $seq =~ /^(\.*)(.*[ACGTNRYMKSWHBVD])(\.*)\s*$/ ) {
      $aStart = length($1);
      $aEnd = $aStart + length($2) - 1;
    }
    my $tSeq = $seq;
    $tSeq =~ s/\.//g;
    my $seqLen = length($tSeq);

    # Mark the sequence as aligned.
    if ( exists $missingSeqIDs{$id} ) {
      delete $missingSeqIDs{$id};
    }
    if ( ! exists $msa{$id} ) {
      $msa{$id} = []; 
    }
    push @{$msa{$id}}, [$start, $end, $aStart, $aEnd, $seqLen, $seq];
  }else {
    warn "Didn't recognize: $_\n";
  }
}
close STK;


#
# Now remove non-colinear fragments, overlaps and collapse individual fragment alignments into one.
# Also creates an insert list for unaligned internal sequences.
#
my @insertList = ();
my %newAlign = ();
foreach my $id ( keys %msa ) {
  if ( @{$msa{$id}} == 1 ){
    $newAlign{$id} = [ $msa{$id}->[0]->[0], $msa{$id}->[0]->[1], $msa{$id}->[0]->[5]];
    next;
  }
  # Sort by length of aligned sequence largest first.
  @{$msa{$id}} = sort { $b->[4] <=> $a->[4] } @{$msa{$id}};
  for ( my $i = 0; $i <= $#{$msa{$id}}; $i++ ){
    my $start = $msa{$id}->[$i]->[0];
    my $end = $msa{$id}->[$i]->[1];
    my $AStart= $msa{$id}->[$i]->[2];
    my $AEnd= $msa{$id}->[$i]->[3];
    for ( my $j = $i+1; $j <= $#{$msa{$id}}; $j++ ){
        my $newStart = $msa{$id}->[$j]->[0];
        my $newEnd= $msa{$id}->[$j]->[1];
        my $newAStart= $msa{$id}->[$j]->[2];
        my $newAEnd= $msa{$id}->[$j]->[3];
        if ( 
             ($newStart >= $start && $newStart <= $end ) || # overlap in sequence space
             ($newEnd >= $start && $newEnd <= $end ) || # overlap in sequence space
             ($newAStart >= $AStart && $newAStart <= $AEnd ) || # overlap in alignment space
             ($newAEnd >= $AStart && $newAEnd <= $AEnd ) || # overlap in alignment space
             ($newAStart > $AStart &&  $newStart < $start )  || # not colinear
             ($newAStart < $AStart && $newStart > $start )      # not colinear
           ) {     # overlap in alignment sapce
          $msa{$id}->[$j]->[0] = "-1";
          $msa{$id}->[$j]->[1] = "-1";
          $msa{$id}->[$j]->[2] = "-1";
          $msa{$id}->[$j]->[3] = "-1";
          $msa{$id}->[$j]->[4] = "-1";
          $msa{$id}->[$j]->[5] = "";
        }
    }
  }
  my $collapsedLine = "";
  my $newStart = 0;
  my $currEnd = 0;
  my $currAEnd = 0;
  @{$msa{$id}} = sort { $a->[0] <=> $b->[0] } @{$msa{$id}};
  for ( my $i = 0; $i <= $#{$msa{$id}}; $i++ ){
    next if ( $msa{$id}->[$i]->[0] == -1 );
    if ( $collapsedLine eq "" ) {
      $newStart = $msa{$id}->[$i]->[0];
      $currEnd = $msa{$id}->[$i]->[1];
      $currAEnd = $msa{$id}->[$i]->[3];
      $collapsedLine = $msa{$id}->[$i]->[5];
      next;
    }
    my $insertionSize = $msa{$id}->[$i]->[0] - $currEnd - 1;
    if ( $insertionSize > 0 ) {
      push @insertList,[$id,$currAEnd+1,substr($seqs{$id},$currEnd,$insertionSize)] 
    }
    $currEnd = $msa{$id}->[$i]->[1];
    $currAEnd = $msa{$id}->[$i]->[3];
    for ( my $j = 0; $j < length($msa{$id}->[$i]->[5]); $j++ ){
      my $newPos = substr($msa{$id}->[$i]->[5],$j,1);
      my $existPos = substr($collapsedLine,$j,1);
      if ( $newPos ne "." ) {
        if ( $existPos ne "." ) {
          warn "Strange...this appears to overlap i=$i!\n";
        }else {
          substr($collapsedLine,$j,1) = $newPos;
        }
      }
    }
  }
  $newAlign{$id} = [ $newStart, $currEnd, $collapsedLine];
}



#print "Dumper: " . Dumper(\@insertList) . "\n";
##print "Dumper: " . Dumper(\%newAlign) . "\n";
#print "<html><pre>\n";
#foreach my $id ( keys %newAlign ) {
#  next if ( $id ne "node-0:100.0" );
#  my $entry = $newAlign{$id};
#  print "" . sprintf("% 15s", $id) . " " . $entry->[2] . "\n";
#}
#exit;


# Put list in reverse order so we can modify sequences without impacting the
# indices
@insertList = sort { $b->[1] <=> $a->[1] } @insertList;
for ( my $i = 0; $i <= $#insertList; $i++ ) {
  my $insID = $insertList[$i]->[0];
  my $alignStart = $insertList[$i]->[1];
  my $seq = $insertList[$i]->[2];
  my $standardLen = 0;
  foreach my $id ( keys %newAlign ) {
    my $entry = $newAlign{$id};
    if ( $id eq $insID ) {
      #substr($entry->[2],$alignStart,0) = "[" . $seq . "]";
      substr($entry->[2],$alignStart,0) = $seq;
    }else {
      #substr($entry->[2],$alignStart,0) = "[" . "."x(length($seq)) . "]";
      substr($entry->[2],$alignStart,0) = "."x(length($seq));
    }
    if ( $standardLen == 0 ) { $standardLen = length($entry->[2]); }
    if ( length($entry->[2]) != $standardLen ) { die "oops...length change at i=$i and alignID=$id\n"; }
    
  }
}

#print "<html><pre style=\"font-family: 'courier new', courier, monospace;\">\n";
#foreach my $id ( keys %newAlign ) {
#  next if ( $id ne "node-0:100.0" );
#  my $s = $newAlign{$id}->[2];
#  $s =~ s/\./-/g;
#  print "" . sprintf("% 15s", $id) . " " . $s . "\n";
#}
#exit;
#print "Dumper: " . Dumper(\%newAlign) . "\n";

my %preSeqs = ();
my %sufSeqs = ();
my $sumPreLens = 0;
my $sumSufLens = 0;
my $coreAlignLen = 0;
foreach my $id ( keys %newAlign ) {
  my $align = $newAlign{$id};
  my $start = $align->[0];
  my $end = $align->[1];
  my $seq = $align->[2];

  if ( $start > 1 ) {
    $preSeqs{$id} = substr($seqs{$id}, 0, $start - 1);
    $sumPreLens += length($preSeqs{$id});
  }
  if ( $end < length($seqs{$id}) ) {
    $sufSeqs{$id} = substr($seqs{$id}, $end);
    $sumSufLens += length($sufSeqs{$id});
  }
  if ( $coreAlignLen == 0 ) {
    $coreAlignLen = length($seq);
  }
}

#
# Add missing sequences in their entirety
#
foreach my $id ( keys(%missingSeqIDs) ) {
  #print STDERR "MISSING: $id\n";
  $preSeqs{$id} = $seqs{$id};
  $sumPreLens += length($seqs{$id});
  $newAlign{$id} = [ -1, -1, "-"x($coreAlignLen)];
}

#print "Sum prefix lengths = $sumPreLens\n";
#print "Sum suffix lengths = $sumSufLens\n";

#print "<html><pre style=\"font-family: 'courier new', courier, monospace;\">\n";

my $nextPrefixIdx = 0;
my $suffixStart = $sumPreLens + $coreAlignLen;
my $nextSuffixIdx = $suffixStart;
foreach my $id ( keys(%newAlign) ){
  my $finalSeq = $newAlign{$id}->[2];
  if ( exists $preSeqs{$id} ) {
    $finalSeq = "-"x($nextPrefixIdx) . $preSeqs{$id} . "-"x($sumPreLens - ($nextPrefixIdx + length($preSeqs{$id}))) . $finalSeq;
    $nextPrefixIdx += length($preSeqs{$id});
  }else {
    $finalSeq = "-"x($sumPreLens) . $finalSeq;
  }
  if ( exists $sufSeqs{$id} ) {
    $finalSeq .= "-"x($nextSuffixIdx - $suffixStart) . $sufSeqs{$id} . "-"x($sumSufLens-(($nextSuffixIdx-$suffixStart)+length($sufSeqs{$id})));
    $nextSuffixIdx += length($sufSeqs{$id});
  }else {
    $finalSeq .= "-"x($sumSufLens);
  }
  $finalSeq =~ s/\./-/g;
  if ( 1 ) {
  print ">$id\n";
  print "$finalSeq\n";
  }else {
    print "" . sprintf("% 15s", $id) . " " . $finalSeq . "\n";
  }
}


