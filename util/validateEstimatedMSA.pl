#!/usr/bin/perl
use strict;

# ./validateEstimatedMSA.pl referenceMSA.fa estimatedMSA.fa or estimatedMSA.stk

open IN,"<$ARGV[0]" or die;
my $id;
my $seq;
my %refSeqs = ();
while ( <IN> ) {
  if ( /^>(\S+)/ ) {
    my $tID = $1;
    if ( $seq ne "" ) {
      $refSeqs{$id} = $seq;
    }
    $id = $tID;
    $seq = "";
    next;
  }
  s/[\n\r\s\-\.]+//g;
  $seq .= uc($_);
}
close IN;
if ( $seq ne "" ) {
  $refSeqs{$id} = $seq;
}


open IN,"<$ARGV[1]" or die;
$id = "";
$seq = "";
my %estSeqs = ();
while ( <IN> ) {
  if ( /^>(\S+)/ ) {
    my $tID = $1;
    if ( $seq ne "" ) {
      $estSeqs{$id} = $seq;
    }
    $id = $tID;
    $seq = "";
    next;
  }

  if ( /^(\S+)\s+([acgtnrymkswhbvdACGTNRYMKSWHBVD\.\-]+)/s*$/ ) {
    my $tID = $1;
    my $tSeq = $2;
    $tSeq =~ s/[\.\-]+//g;
    $estSeqs{$tID} = uc($tSeq);
    $id = "";
    $seq = "";
    next;
  }
 
  if ( /^([acgtnrymkswhbvdACGTNRYMKSWHBVD\.\-]+)\s*$/ ) {
    my $tSeq = $1;
    $tSeq =~ s/[\.\-]+//g;
    $seq .= uc($tSeq);
  }
}
if ( $seq ne "" ) {
  $estSeqs{$id} = $seq;
}
close IN;


my $errors = 0;
foreach my $refID ( keys(%refSeqs) ){
  if ( exists $estSeqs{$refID} ) {
    if ( $estSeqs{$refID} ne $refSeqs{$refID} ) {
      my $firstDiffPos = -1;
      for ( my $i = 0; $i < length($refSeqs{$refID}); $i++ ){
        if ( $i == length($estSeqs{$refID}) ) {
          $firstDiffPos = $i;
          last;
        }
        if ( substr($estSeqs{$refID},$i,1) ne substr($refSeqs{$refID},$i,1) ) {
           $firstDiffPos = $i; 
           last;
        }
      }
      print "Sequence $refID differs: refLen = " . length($refSeqs{$refID}) . " estLen = " .  length($estSeqs{$refID}) . " firstDiff @ " . $firstDiffPos . ":\n";
      print "    ref: $refSeqs{$refID}\n";
      print "    est: $estSeqs{$refID}\n";
      $errors++;
    }
  }else {
    print "Sequence $refID is missing from estimated alignment!\n";
    $errors++;
  }
}

if ( $errors ne 0 ) {
  die "There are $errors errors in $ARGV[1]!\n";
}

