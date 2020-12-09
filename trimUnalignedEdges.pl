#!/usr/local/bin/perl
use strict;


open IN,"$ARGV[0]" or die;
my @ids = ();
my @seqs = ();
my $id;
my $seq;
while ( <IN> ) {
  if ( /^>(\S+)/ ) {
    my $tmpID = $1;
    if ( $seq ) {
      push @seqs, $seq;
      push @ids, $id;
    }
    $id = $tmpID;
    $seq = "";
    next;
  }
  s/[\n\r\s]+//g;
  $seq .= $_;
}
if ( $seq ) {
  push @seqs, $seq;
  push @ids, $id;
}
close IN;


my $start = -1;
for ( my $i = 0; $i < length($seqs[0]); $i++ ) {
  my $coverage = 0;
  for ( my $j = 0; $j <= $#seqs; $j++ ) {
    if ( substr($seqs[$j],$i,1) ne "-" ) {
      $coverage++;
    }
  }
  if ( $start < 0 && $coverage > 1 ){
    $start = $i;
    last;
  }
}

my $end = -1;
for ( my $i = length($seqs[0])-1; $i >= 0; $i-- ) {
  my $coverage = 0;
  for ( my $j = 0; $j <= $#seqs; $j++ ) {
    if ( substr($seqs[$j],$i,1) ne "-" ) {
      $coverage++;
    }
  }
  if ( $end < 0 && $coverage > 1 ) {
    $end = $i;
    last;
  }
}

#print "Start = $start and end = $end\n";
for ( my $i = 0; $i <= $#seqs; $i++ ) {
  my $seq = substr($seqs[$i], $start, $end - $start + 1);
  if ( $seq =~ /^-+$/ ) {
  }else {
  print ">$ids[$i]\n";
  print "" . substr($seqs[$i], $start, $end - $start + 1) . "\n";
  }
}


