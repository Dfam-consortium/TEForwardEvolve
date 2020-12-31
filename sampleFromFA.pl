#!/usr/local/bin/perl
use strict;
use File::Basename;

my $sampleSizes= $ARGV[0];
my $fastaFile = $ARGV[1];
if ( @ARGV > 2 ) {
  srand($ARGV[2]);
}

my ($faName, $faPath, $faSuffix) = fileparse($fastaFile);
$faName = $1 if ( $faName =~ /^(\S+)\..*/ );

my $id;
my $seq;
my @seqs;
open IN,"<$fastaFile" or die;
while (<IN>) {
  if ( /^>(\S+)/ ) {
    my $tID = $1;
    if ( $seq && length($seq) >= 50 ) {
      push @seqs,[$id,$seq];
    }
    $id = $tID;
    $seq = "";
    next;
  }
  s/[\n\r\s]+//g;
  $seq .= $_;
}
if ( $seq && length($seq) >= 50 ) {
  push @seqs,[$id,$seq];
}
close IN;

$sampleSizes =~ s/[\s\"\']//g;
my @sizes = split(/,/,$sampleSizes);
foreach my $size ( @sizes ) {
  last if ( $size > scalar @seqs );
  my $outFile = "$faPath/$faName-$size"."sample.fa";
  open OUT,">$outFile" or die;
  my $numColl = 0;
  my %dup = ();
  while ( 1 ) {
    my $idx = int(rand($#seqs));
    if ( ! exists $dup{$idx} ) {
      $dup{$idx}++;
      #print OUT ">$seqs[$idx]->[0]\n$seqs[$idx]->[1]\n";
      print OUT ">seq-$idx  $seqs[$idx]->[0]\n$seqs[$idx]->[1]\n";
      $numColl++;
      last if ( $numColl >= $size );
    }
  }
  close OUT;
}


1;
