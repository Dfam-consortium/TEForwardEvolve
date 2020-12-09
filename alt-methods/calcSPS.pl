#!/usr/local/bin/perl
#
# calcSPS.pl : Very slow perl implementation of the SPS calculation
#
# -Robert Hubley 2020
#
use strict;
use Data::Dumper;


my $testSeqs = readFASTA($ARGV[0]);
my $refSeqs = readFASTA($ARGV[1]);
#my $testSeqs = readFASTA("refiner.fa");
#my $testSeqs = readFASTA("muscle.fa");


my %refColToBase;
my %testBaseToCol;
my $numRefCols = 0;
foreach my $rID ( keys(%{$refSeqs}) ) {
  my $c_rSeq = $refSeqs->{$rID};
  $numRefCols = length($c_rSeq) if ( $numRefCols == 0 );
  $c_rSeq =~ s/-//g;
 
  my $c_tSeq = $testSeqs->{$rID};
  $c_tSeq =~ s/-//g;

  my $start = index($c_rSeq, $c_tSeq);

  if ( $start < 0 ) {
    print "Ooops...rID=$rID c_rSeq=$c_rSeq ctSeq=$c_tSeq\n";
  }
 
  # Handle the fact that a multiple alignment may not align all edge bases
  for ( my $i = 0; $i < $start; $i++ ){
    $testBaseToCol{$rID}->[$i] = -1;
  }
  my $residueIdx = $start;
  for ( my $i = 0; $i < length($testSeqs->{$rID}); $i++ )
  {
    my $chr = substr($testSeqs->{$rID},$i,1);
    if ( $chr !~ /[\.\- ]/ ) {
      $testBaseToCol{$rID}->[$residueIdx] = $i;
      $residueIdx++;
    }
  }
  for ( my $i = $residueIdx; $i < length($c_tSeq); $i++ ){
    $testBaseToCol{$rID}->[$i] = -1;
  }
  
  $residueIdx = 0;
  for ( my $i = 0; $i < length($refSeqs->{$rID}); $i++ ) 
  {
     my $chr = substr($refSeqs->{$rID},$i,1);
     if ( $chr !~ /[\.\- ]/ ) {
       $refColToBase{$rID}->[$i] = $residueIdx;
       $residueIdx++;
     }else {
       $refColToBase{$rID}->[$i] = -1;
     }
  }
}

#print "Dumper : " . Dumper(\%refColToBase) . "\n";
#print "Dumper : " . Dumper(\%testBaseToCol) . "\n";
   

my $DEBUG = 0;
my $SCScore = 0;
my $SRScore = 0;
for ( my $i = 0; $i < $numRefCols; $i++ ){
  my $colSum = 0;
  my $refSum = 0;
  foreach my $rID1 ( keys(%{$refSeqs}) ) {
    if ( $refColToBase{$rID1}->[$i] >= 0 )
    {
      # This is occupied with base #
      my $refBaseBasePos1 = $refColToBase{$rID1}->[$i];
      my $equivTestBaseCol1 = $testBaseToCol{$rID1}->[$refBaseBasePos1];
      foreach my $rID2 ( keys(%{$refSeqs}) ) {
        if ( $rID1 ne $rID2 && $refColToBase{$rID2}->[$i] >= 0 )
        {
print "$i: $rID1, $rID2\n";
          $refSum++;
          # This is occupied with base #
          my $refBaseBasePos2 = $refColToBase{$rID2}->[$i];
          my $equivTestBaseCol2 = $testBaseToCol{$rID2}->[$refBaseBasePos2];
          $colSum += 1 if ( $equivTestBaseCol1 >= 0 && $equivTestBaseCol1 == $equivTestBaseCol2 );
        }
      }
    }
  }
  $SCScore += $colSum;
  $SRScore += $refSum;
} # for ( my $i
print "SPS Score = $SCScore / $SRScore = " . ( $SCScore / $SRScore ) . "\n";


sub readFASTA {
  my $file = shift;

  my %seqs;

  my $id; 
  my $seq;
  open IN,"<$file" or die "Could not open $file for reading!\n";
  while( <IN> ) {
    if (/^>(\S+)/) {
    my $tID = $1;
      if ($seq) {
        $seqs{$id} = uc($seq);
      }
      $id = $tID;
      $seq = "";
      next;
    } 
    s/[\n\r\s]//g;
    $seq .= $_;
  }
  if ( $seq ) {
    $seqs{$id} = uc($seq);
  }
  return \%seqs;
}


