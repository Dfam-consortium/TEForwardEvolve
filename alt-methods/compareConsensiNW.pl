#!/usr/local/bin/perl
use strict;
use lib "/usr/local/RepeatMasker";
use lib "/home/rhubley/projects/RepeatModeler";
use NeedlemanWunschGotohAlgorithm;
use SequenceSimilarityMatrix;
use SearchResult;

my $lupMatrix = SequenceSimilarityMatrix->new();
$lupMatrix->parseFromFile( "/usr/local/RepeatModeler/matrices/ncbi/nt/comparison.matrix" );

my $predictedFA = $ARGV[0];
my $referenceFA = $ARGV[1];

print "compareConsensiNW.pl\n";
print "Predicted Consensus: $predictedFA\n";
print "Reference Consensus: $referenceFA\n";

my ($pred_ID, $pred ) = readSingleSeqFasta( $predictedFA );
my ($ref_ID, $ref ) = readSingleSeqFasta( $referenceFA );

my $sr = NeedlemanWunschGotohAlgorithm::search(
  querySeq   => uc($pred),
  subjectSeq => uc($ref),
  matrix         => $lupMatrix,
  insOpenPenalty => -25,
  insExtPenalty  => -5,
  delOpenPenalty => -25,
  delExtPenalty  => -5
);

$sr->setQueryName("predicted");
$sr->setSubjName("reference");
$sr->setMatrixName("comparison.matrix");

my $totalRefCpGs = 0;
my $restoredCpGs = 0;
my $predAln = $sr->getQueryString();
my $refAln = $sr->getSubjString();
for( my $i = 0; $i < length($refAln); $i++ ){
  my $rbase1 = substr($refAln,$i,1);
  if ( $rbase1 eq "C" ) {
    my $rbase2 = "";
    for( my $j = $i+1; $j < length($refAln); $j++ ){
      $rbase2 = substr($refAln,$j,1);
      next if ( $rbase2 eq "-" );
      if ( $rbase2 eq "G" ) {
        $totalRefCpGs++;
        if ( substr($predAln,$i,1) eq "C" && substr($predAln,$j,1) eq "G" ) {
          $restoredCpGs++;
        }
      }
      last;
    }
  }
}

my $rawScore = $sr->getScore();
my $matrix =
      Matrix->new(
        fileName => "/usr/local/RepeatModeler/matrices/ncbi/nt/comparison.matrix"  );


my ( $rescore, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions ) = $sr->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -25,
                             gapExtPenalty => -5,
                             divCpGMod => 1,
                             complexityAdjust => 1
                             );

$sr->setScore($rescore);
$sr->setPctDiverge( sprintf( "%0.2f", $divergence ) );
$sr->setPctInsert( sprintf( "%0.2f",  $percIns ) );
$sr->setPctDelete( sprintf( "%0.2f",  $percDel ) );

my $predLen = length($pred);
my $tpred = $pred;
$tpred =~ s/N//g;
my $predLenNonAmbig = length($tpred);

my $predAln = $sr->getQueryString();
my $refAln = $sr->getSubjString();

my ($delCnt) = ( $predAln =~ tr/-/-/ );
my ($insCnt) = ( $refAln =~ tr/-/-/ );


print "" . $sr->toStringFormatted(SearchResult::AlignWithQuerySeq) . "\n";
print "Raw Score = $rawScore\n";
print "Complexity Adjusted Score = $rescore\n";
print "Kimura = " . sprintf("%0.2f",$divergence) . " %\n";
print "Percent Substitutions = " .  sprintf("%0.2f",(100 *(( $transitions + $transversions ) / $well_characterized_bases))) . " %\n";
print "  - Transitions = $transitions (CpG adj)\n";
print "  - Transversions = $transversions\n";
print "Percent Insertion = " . sprintf("%0.2f",($insCnt/length($ref))*100) . " % [ $insCnt bp ]\n";
print "Percent Deletion = " . sprintf("%0.2f",($delCnt/length($ref))*100) . " % [ $delCnt bp ]\n";
if ( $totalRefCpGs > 0 ) {
print "Reference CpG sites predicted: $restoredCpGs/$totalRefCpGs ( " .  sprintf("%0.2f",($restoredCpGs/$totalRefCpGs)*100) . " % )\n";
}else {
print "Reference CpG sites predicted: 0/0 ( 100 % )\n";
}
print "predicted cons length = $predLen [ $predLenNonAmbig non-ambiguous bp ]\n";
print "reference cons length = " . length($ref) . "\n";


sub readSingleSeqFasta {
  my $file = shift;

  my $id;
  my $seq;
  open IN,"<$file" or die "could not open $file\n";
  while ( <IN> ) {
    if ( /^>(\S+)/ ) {
      my $tmpID = $1;
      if ( $seq ) {
        last;
      }
      $id = $tmpID;
      $seq = "";
      next;
    }
    s/[\n\r\s]+//g;
    $seq .= $_;
  }
  close IN;
  return($id, $seq);
}

