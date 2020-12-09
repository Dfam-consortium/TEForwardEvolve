#!/usr/bin/perl
use strict;

my $scaleFactor = 0.000001;


my @matrix = ();
my $parse = 0;
while ( <> ) {
  # RATE_MAT:d
  if ( /RATE_MAT:/ ) {
    $parse = 1;
    next;
  }
  if ( /^TREE:/ ) { last; }
  if ( $parse && /^\s*([\d\-\.]+.*)/ ) {
    my @flds = split(/\s+/,$1);
    if ( @flds == 64 ) {
      push @matrix,[@flds];
    }else {
      warn "Doesn't look like matrix data: $_\n";
    }
  }
}

my %tripletHash = ();
my @tripletArray = ();
my $index = 0;
foreach my $b1 ( 'A','C','G','T' ) {
  foreach my $b2 ( 'A','C','G','T' ) {
    foreach my $b3 ( 'A','C','G','T' ) {
      my $triplet = $b1.$b2.$b3;
      push @tripletArray,$triplet;
      $tripletHash{$triplet} = $index;
      $index++;
    }
  }
}

#trinucleotide   A       C       G       T
print "trinucleotide\tA\tC\tG\tT\n";
$index = 0;
foreach my $triplet ( @tripletArray ) {
  my $mutTri = $triplet;
  print "$triplet";
  foreach my $sub ( 'A','C','G','T' ) {
    substr($mutTri,1,1) = $sub;
    if ( $mutTri eq $triplet ) {
      print "\t0";
    }else {
      print "\t" . (($matrix[$index]->[$tripletHash{$mutTri}])*$scaleFactor);
    }
  }
  print "\n";
  $index++;
}


    
