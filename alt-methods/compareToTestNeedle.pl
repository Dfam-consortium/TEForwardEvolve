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
my $testSeqsFA = $ARGV[1];

print "compareToTestNeedle.pl\n";
print "Predicted Consensus: $predictedFA\n";
print "Test Set Sequences: $testSeqsFA\n";

if ( ! -s $predictedFA ) {
  print "Cannot score consensi: predicted fasta file $predictedFA is empty!\n";
  exit;
}
if ( ! -s $testSeqsFA ) {
  print "Cannot score consensi: test set fasta file $testSeqsFA is empty!\n";
  exit;
}


open IN,"$needle -datafile $matrixFile -gapopen $gapOpen -gapextend $gapExtn -endweight 1 -endopen $gapOpen -endextend $gapExtn -aformat3 score -outfile stdout $predictedFA $testSeqsFA |" or die "Could not execute needle!\n";
my $score;
my $cnt = 0;
while (<IN>) {
  if ( /^(\S+)\s+(\S+)\s+(\d+)\s+\(([\-\d\.]+)\)/ ) {
    $score += $4;
    $cnt++;
  }
}
close IN;

print "Alignments = $cnt\n";
print "Sum Score = $score\n";
print "Average Score = " . ($score/$cnt) . "\n";
