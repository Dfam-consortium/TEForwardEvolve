#!/usr/local/bin/perl
use strict;
use lib "/usr/local/RepeatMasker";
use SearchResult;
use Matrix;

my $EMBOSSDIR = "/usr/local/EMBOSS-6.6.0/bin";
my $needle = "$EMBOSSDIR/needle";
#my $matrixFile = "/usr/local/RepeatModeler/matrices/ncbi/nt/comparison.matrix";
my $matrixFile = "/home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/comparison-deterNs.matrix";
my $gapOpen = 25.0;
my $gapExtn = 5.0;

my $predictedFA = $ARGV[0];
my $referenceFA = $ARGV[1];

print "compareConsensiNeedle.pl\n";
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


open IN,"$needle -datafile $matrixFile -gapopen $gapOpen -gapextend $gapExtn -endweight 1 -endopen $gapOpen -endextend $gapExtn -outfile stdout $predictedFA $referenceFA |" or die "Could not execute needle!\n";
my $score;
my @seqs = ();
my $flip = 0;
while (<IN>) {
  if ( /^#\s+Score:\s+([\-\d\.]+)/ ) {
    $score = $1;
  }
  if ( /^(\S+)\s+\d+\s+([ACGTRYSWKMBDHVNacgtryswkmbdhvn\-]+)\s+\d+\s*$/ ) {
    $seqs[$flip] .= uc($2);
    $flip ^= 1;
  }
}
close IN;

my $tStr = $seqs[0];
$tStr =~ s/-//g;
my $qLen = length($tStr);
$tStr = $seqs[1];
$tStr =~ s/-//g;
my $sLen = length($tStr);

my $rawScore = $score;

my $sr = SearchResult->new( score => $score,
                            queryName => "predicted",
                            subjName => "reference",
                            queryStart => 1,
                            queryEnd => $qLen,
                            queryString => $seqs[0], 
                            orientation => "",
                            subjStart => 1,
                            subjEnd => $sLen,
                            subjString => $seqs[1]);

$sr->setMatrixName("comparison.matrix");


my $pred = $seqs[0];
my $ref = $seqs[1];
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

$pred =~ s/-//g;
$ref =~ s/-//g;
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


