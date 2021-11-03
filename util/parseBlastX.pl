#!/usr/local/bin/perl
use strict;

#print "Curated Proteins\n";
#print "================\n";
foreach my $SIM ("Arthur2", "Tigger10", "Zaphod", "Zaphod2") {
  opendir DIR,"$SIM-eval" or die;
  while ( my $entry = readdir(DIR) ) {
    # Curated proteins
    # orffrags-500sample-mafft.blastx
    # Homologs
    # orffrags-500sample-mafft.dbprot.blastx
    if ( $entry =~ /orffrags-150sample-(\S+)\.dbprot.blastx/ ) {
      #doit("$SIM-eval/$entry",$SIM,$1);
    }elsif ( $entry =~ /orffrags-150sample-(\S+)\.blastx/ ) {
      doit("$SIM-eval/$entry",$SIM,$1);
    }
  }
}
#print"\n\n";
#print "Homolog Proteins\n";
#print "================\n";
foreach my $SIM ("Arthur2", "Tigger10", "Zaphod", "Zaphod2") {
  opendir DIR,"$SIM-eval" or die;
  while ( my $entry = readdir(DIR) ) {
    # Curated proteins
    # orffrags-500sample-mafft.blastx
    # Homologs
    # orffrags-500sample-mafft.dbprot.blastx
    if ( $entry =~ /orffrags-150sample-(\S+)\.dbprot.blastx/ ) {
      doit("$SIM-eval/$entry",$SIM,$1);
    }elsif ( $entry =~ /orffrags-150sample-(\S+)\.blastx/ ) {
      #doit("$SIM-eval/$entry",$SIM,$1);
    }
  }
}


sub doit {
  my $file = shift;
  my $family = shift;
  my $method = shift;
my $proteinID = "";
my $prevID = "";
my $expect = "";
my $prevExpect = "";
my $protLen = 0;
my $prevProtLen = 0;
my $prevLine = "";
my $inAlign = 0;
my @ranges = ();
open IN,"<$file" or die;
while ( <IN> ) {
  if ( /^>\s*(\S+)/ ) {
    $proteinID = $1;
  }

  if ( /^Length=(\d+)/ ) {
    $protLen = $1;
  }

  if ( /,\s+Expect = (\S+),/ ||
       /,\s+Expect\(\d+\)\s+=\s+(\S+),/ ) {
    $expect = $1;
    if ( $prevExpect ne "" && $expect ne $prevExpect ) {
      if ( @ranges ) {
      #if ( @ranges  && $prevExpect < 0.001 ) {
        print "$family\t$prevProtLen\t$method\t$prevID\t$prevExpect";
        foreach my $range ( @ranges ) {
          print "\t[" . $range->[0] . "," . $range->[1] . "]";
        }
        print "\n";
      }
      @ranges = ();
    }
    $prevID = $proteinID;
    $prevExpect = $expect;
    $prevProtLen = $protLen;
 }

  if ( /Sbjct\s+(\d+)\s+\S+\s+(\d+)/ ) 
  {
    if ( @ranges && abs($ranges[-1]->[1] - $1) == 1 )
    {
      $ranges[-1]->[1] = $2;
    }else {
      push @ranges, [$1,$2];
    }
  }
}
close IN;

    #if ( @ranges  && $prevExpect < 0.001 ) {
    if ( @ranges ) {
      print "$family\t$prevProtLen\t$method\t$prevID\t$prevExpect";
      foreach my $range ( @ranges ) {
        print "\t[" . $range->[0] . "," . $range->[1] . "]";
      }
      print "\n";
    }
}
 
#> XP_012561266.1 PREDICTED: zinc finger MYM-type protein 1-like 
#[Hydra vulgaris]
#Length=819
#
# Score = 39.3 bits (90),  Expect(2) = 1e-07, Method: Compositional matrix adjust.
# Identities = 26/62 (42%), Positives = 36/62 (58%), Gaps = 5/62 (8%)
# Frame = +2
#
#Query  344  LFQRPKFNQLEXFXKXHPIQPFET-KNLPF-NGKSAFHR-KDGTXCPWLSY--SPEXVFY  508
#            +F RP+   L+ F K HP+QP E  +NLPF N +  F+R +D     WLSY  S + +F 
#Sbjct  98   IFMRPQSQCLKSFFKSHPVQPKEDGENLPFSNVRRVFYRSQDDINRVWLSYCHSSKAMFC  157
#
#Query  509  RV  514
#             V
#Sbjct  158  TV  159
#
#
# Score = 22.7 bits (47),  Expect(2) = 1e-07, Method: Compositional matrix adjust.
# Identities = 10/32 (31%), Positives = 14/32 (44%), Gaps = 0/32 (0%)
# Frame = +1
#
#Query  490  SRXGFLPXASAYSXDLDFSKFISDXKDWRHSH  585
#            S+  F     AYS   + + F     DW+H H
#Sbjct  152  SKAMFCTVCLAYSKASESNAFTDGMNDWKHVH  183
