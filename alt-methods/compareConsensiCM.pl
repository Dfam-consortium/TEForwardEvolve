#!/usr/local/bin/perl
use strict;
use lib "/usr/local/RepeatMasker";
use SearchResult;
use Matrix;

my $FASTADIR = "/usr/local/fasta";
my $ggsearch = "$FASTADIR/ggsearch36";
my $matrixFile = "/usr/local/RepeatModeler/matrices/ncbi/nt/comparison.matrix";
my $gapOpen = -25;
my $gapExtn = -5;

my $predictedFA = $ARGV[0];
my $referenceFA = $ARGV[1];

print "compareConsensi.pl\n";
print "Predicted Consensus: $predictedFA\n";
print "Reference Consensus: $referenceFA\n";

if ( ! -s $predictedFA ) {
  print "Cannot score consensi: predicted fasta file $predictedFA is empty!\n";
  exit;
}
if ( ! -s $referenceFA ) {
  print "Cannot score consensi: reference fasta file $referenceFA is empty!\n";
  exit;
}




open IN,"$ggsearch -n -s $matrixFile -f $gapOpen -g $gapExtn -m 10 $referenceFA $predictedFA|" or die;
my $inAlign = 0;
my $start = -1;
my $stop = -1;
my $alignData = "";
my $score;
my @seqs = ();
my $prevScore = 0;
while (<IN>) {
  $score = $1 if ( /^;\s+gnw_score:\s+([\d\-]+)/ );
  $start = $1 if ( /^;\s+al_start:\s+(\d+)/ );
  $stop = $1 if ( /^;\s+al_stop:\s+(\d+)/ );
  if ( /^;\sal_display_start/ ) {
    $inAlign = 1;
    next;
  }
  if ( $inAlign && /^[>;]/ ) {
    if ( $score > $prevScore ) {
      $prevScore = $score;
      @seqs = ();
    }
    push @seqs,[$start, $stop, $alignData, $score];
    $alignData = "";
    $inAlign = 0;
    $start = -1;
    $stop = -1;
    next;
  }
  if ( $inAlign ) {
    s/[\n\r\s]+//g;
    $alignData .= $_;
    next;
  }
}
close IN;

if ( @seqs ) {
my $orient = "";
$orient = "C" if ( $seqs[0]->[0] > $seqs[0]->[1] );

my $sr;
if ( $orient eq "C" ) {
 $sr = SearchResult->new( score => $seqs[0]->[3],
                            queryName => "predicted",
                            subjName => "reference",
                            queryStart => $seqs[1]->[0],
                            queryEnd => $seqs[1]->[1],
                            queryString => $seqs[1]->[2], 
                            orientation => $orient,
                            subjStart => $seqs[0]->[1],
                            subjEnd => $seqs[0]->[0],
                            subjString => $seqs[0]->[2]);
}else {
 $sr = SearchResult->new( score => $seqs[0]->[3],
                            queryName => "predicted",
                            subjName => "reference",
                            queryStart => $seqs[1]->[0],
                            queryEnd => $seqs[1]->[1],
                            queryString => $seqs[1]->[2], 
                            orientation => $orient,
                            subjStart => $seqs[0]->[0],
                            subjEnd => $seqs[0]->[1],
                            subjString => $seqs[0]->[2]);
}


my $pred = $seqs[1]->[2];
my $ref = $seqs[0]->[2];
my $totalRefCpGs = 0;
my $restoredCpGs = 0;
for( my $i = 0; $i < length($ref); $i++ ){
  my $rbase1 = substr($ref,$i,1);
  if ( $rbase1 eq "C" ) {
    my $rbase2 = "";
    for( my $j = $i+1; $j < length($ref); $j++ ){
      $rbase2 = substr($ref,$j,1);
      next if ( $rbase2 eq "-" );
      if ( $rbase2 eq "G" ) {
        $totalRefCpGs++;
        if ( substr($pred,$i,1) eq "C" && substr($pred,$j,1) eq "G" ) {
          $restoredCpGs++;  
        }
      }
      last;
    }
  }
}

my $matrix =
      Matrix->new(
        fileName => $matrixFile );


my ( $rescore, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions ) = $sr->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -25,
                             gapExtPenalty => -5,
                             divCpGMod => 1,
                             complexityAdjust => 1 );

$sr->setScore($rescore);
$sr->setPctDiverge( sprintf( "%0.2f", $divergence ) );
$sr->setPctInsert( sprintf( "%0.2f",  $percIns ) );
$sr->setPctDelete( sprintf( "%0.2f",  $percDel ) );

print "" . $sr->toStringFormatted(SearchResult::AlignWithQuerySeq) . "\n";
print "Score = $rescore\n";
print "Kimura = " . sprintf("%0.2f",$divergence) . " %\n";
print "Percent Substitutions = " .  sprintf("%0.2f",(100 *(( $transitions + $transversions ) / $seqs[0]->[1]))) . " %\n";
print "  - Transitions = $transitions\n";
print "  - Transversions = $transversions\n";
print "Percent Insertion = " . sprintf("%0.2f",$percIns) . " %\n";
print "Percent Deletion = " . sprintf("%0.2f",$percDel) . " %\n";
if ( $totalRefCpGs > 0 ) {
print "Reference CpG sites predicted: $restoredCpGs/$totalRefCpGs ( " .  sprintf("%0.2f",($restoredCpGs/$totalRefCpGs)*100) . " % )\n";
}else {
print "Reference CpG sites predicted: 0/0 ( 100 % )\n";
}
}else {
  print "Cannot score consensi: No alignment returned by ggsearch\n";
}
                            
