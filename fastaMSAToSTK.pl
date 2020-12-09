#!/usr/local/bin/perl
use strict;

open FA,"<$ARGV[0]" or die;


print "# STOCKHOLM 1.0\n";
print "#=GF ID $ARGV[0]\n";
my $seq;
my $id;
while (<FA>) {
  if ( />(\S+)/ ) {
    my $tmpID = $1;
    if ( $tmpID =~ /.*(node-\d+:\d+\.\d+).*/ ) {
      $tmpID = $1;
    }
    if ( $seq ) {
      print "$id $seq\n";
    }
    $id = $tmpID;
    $seq = "";
    next;
  }
  s/[\n\r\s]+//g;
  s/\-/\./g;
  $seq .= $_;
}
if ( $seq ) {
  print "$id $seq\n";
}
print "//\n";

#cat mus.stk | perl -ne '{ if ( /^Unknown:(node-\d+:\d+\.\d+):\S+\s+(\S+)/ ) { print "$1 $2\n"; }else { print; }}' > mus2.stk

#cat ../../TEForwardEvolve/LINETree-1/rep-1/gput4000-arian.fa | perl -ne '{ if ( $start ne "1" ) { print "# STOCKHOLM 1.0\n"; print "#=GF ID foo\n"; $start = "1" } if ( />(\S+)/ ) { $seqID = $1; next; } s/\-/\./g; s/[\s\n\r]+//g; $seq = $_;  print "$seqID $seq\n";  if ( eof() ) { print "//\n"; }}' > arian-full.stk








