#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) TEForwardEvolve.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A DNA sequence evolver with a focus on Transposable Element
##      sequences. 
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
#
=head1 NAME

TEForwardEvolve.pl - A DNA sequence evolver with triplet rate matrices and indels

=head1 SYNOPSIS

  TEForwardEvolve.pl [-version] [-tree <newick_file>] [-seed <fasta_file>] 
                     [-matrix <matrix_file>]
                     [-generations_per_unit_time <factor>] 
                     [-insertion_rate <rate>][-deletion_rate <rate>]
                     [-mean_del_size <bp>][-mean_ins_size <bp>]
                     [-max_del_size <bp>][-max_ins_size <bp>]
                     [-min_frag_len <bp>][-verbosity #]
                     [-srand #]

=head1 DESCRIPTION

This sequence evolver combines two aspects of sequence evolution not found
currently in the same tree-based forward evolver:

  - Triplet rate matrices to more accurately model flanking base effects
    on substitution including CpG mutation rates in mammalian DNA.

  - An indel model.

INDELible and indel-seq-gen both provide the combination of substition and indel
models but neither support triplet matrices.  Trevolver supports triplet matrices
but not an indel model.

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2020 Robert Hubley, Institute for Systems Biology

=head1 LICENSE

This code may be used in accordance with the Open Source License v. 3.0
http://opensource.org/licenses/OSL-3.0

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin;
use POSIX;

my $Version = "0.1";

# Globals for BFS algorithms
my %globalNodeNamesHash = ();
my @globalCombinedMuts = ();
my @globalNodeList = ();
my $verbosity = 0;
# Globals for DFS_avgSubs
my $extantSubSum = 0;
my $extantCnt = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-version', # print out the version and exit
    '-tree=s',
    '-seed=s',
    '-matrix=s',
    '-srand=s',
    '-mutation_rate=s',
    '-transition_factor=s',
    '-cpg_factor=s',
    '-generations_per_unit_time=s',
    '-insertion_rate=s',
    '-deletion_rate=s',
    '-mean_del_size=s',
    '-max_del_size=s',
    '-mean_ins_size=s',
    '-max_ins_size=s',
    '-verbosity=s'
);

my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions(\%options, @getopt_args)) {
    usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ($options{'version'}) {
  print "$Version\n";
  exit;
}

my $newickFile = "tree.nw";
if ( exists $options{'tree'} && -s $options{'tree'} ) {
  $newickFile = $options{'tree'};
}

# Read Newick tree
#
#   S3
#   /\
# s1 s2
#
# (s1:#,s2:#)s3;
# 
open IN,"<$newickFile" or die "ERROR: Could not open Newick file $newickFile!\n";
my $nw_data;
while ( <IN> ) {
  s/[\n\r\s]+//g;
  $nw_data .= $_;
}
close IN;

# Based loosely on a little format conversion trick by
# Kyle Hasselbacher ( https://www.perlmonks.org/?node_id=717769 ).
# Modified to work with named subnodes, and non-binary structures.
#
# Generates records like { 'VALUE' => s1:0.1, 
#                          'CHILDREN' => [ {'VALUE' => "s2:0.1", ... }, 
#                                          {'VALUE' => "s3:0.1", ... } ]
#                        }
$nw_data =~ s/([^,()]+)/\{VALUE=>"$1"\}/g;
$nw_data =~ s/\(/\{CHILDREN=>\[/g;
$nw_data =~ s/\){(VALUE[^}]+)}/\],$1}/g;
$nw_data =~ s/\)/\}/g;
my $tref = eval $nw_data;
undef $nw_data;

#
# Fix Missing Node Names and initialize a global
# list of nodes.
#
breadthFirstSearch($tref, \&fixMissingNodeNames);


# 
# Read in the seed sequence
#
my $seedFile = "seed.fa";
if ( exists $options{'seed'} && -s $options{'seed'} ) {
  $seedFile = $options{'seed'};
}
open IN,"$seedFile" or die "ERROR: Could not open $seedFile for reading!\n";
my $seedSeq;
while ( <IN> ) {
  if ( /^>/ ) {
    if ( $seedSeq ) {
      last;
    }
    next;
  } 
  s/[\n\r\s]+//g;
  if ( /[^ACGTacgt]/ ) {
    die "ERROR: Seed sequence contain ambiguous base codes.  Must be only A,C,G,T! '$_'\n";
  }
  $seedSeq .= uc($_);
}
close IN;


#
# Setup initial root sequence
#
$tref->{'SEQUENCE'} = $seedSeq;

my @triAlphabet = (
'AAA', 'AAC', 'AAG', 'AAT', 'ACA',
'ACC', 'ACG', 'ACT', 'AGA', 'AGC',
'AGG', 'AGT', 'ATA', 'ATC', 'ATG',
'ATT', 'CAA', 'CAC', 'CAG', 'CAT',
'CCA', 'CCC', 'CCG', 'CCT', 'CGA',
'CGC', 'CGG', 'CGT', 'CTA', 'CTC',
'CTG', 'CTT', 'GAA', 'GAC', 'GAG',
'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
'GGA', 'GGC', 'GGG', 'GGT', 'GTA',
'GTC', 'GTG', 'GTT', 'TAA', 'TAC',
'TAG', 'TAT', 'TCA', 'TCC', 'TCG',
'TCT', 'TGA', 'TGC', 'TGG', 'TGT',
'TTA', 'TTC', 'TTG', 'TTT' );

# Parameters for theoretical matrices
my $cpgFactor = 20;
$cpgFactor = $options{'cpg_factor'} if ($options{'cpg_factor'} );

my $transIFactor = 1;
$transIFactor = $options{'transition_factor'} if ($options{'transition_factor'} );

my $mutRate = 7.5E-07;
$mutRate = $options{'mutation_rate'} if ($options{'mutation_rate'} );

#
# Read in substitution matrix
#
my $triMatrix;
if ( $options{'matrix'} ) {
  # Custom Matrix
  $triMatrix = {};
  open IN,"<$options{'matrix'}" or die "Could not open matrix file $options{'matrix'}\n";
  while ( <IN> ) {
     #trinucleotide   A       C       G       T
     next if ( /^\s*trinucleotide\s+A\s+C\s+G\s+T\s*$/ );
     next if ( /^\s*$/ );
     if ( /^([ACGTacgt]{3,3})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/ ) {
       $triMatrix->{$1}->{'A'} = $2;
       $triMatrix->{$1}->{'C'} = $3;
       $triMatrix->{$1}->{'G'} = $4;
       $triMatrix->{$1}->{'T'} = $5;
     }else {
       die "Error in matrix file.  Didn't expect:\n$_\n\nMatrix file should look like:\n" .
           "trinucleotide   A       C       G       T\n" .
           "AAA     0       2.50E-07        2.50E-07        2.50E-07\n" .
           "AAC     0       2.50E-07        2.50E-07        2.50E-07\n" .
           "AAG     0       2.50E-07        2.50E-07        2.50E-07\n" .
           "...\n\n";
     }
  }
  close IN;
  # Validate matrix
  foreach my $triWord ( @triAlphabet ){
    if ( not exists $triMatrix->{$triWord} ) {
      die "Matrix is missing tri-nucleotide word: $triWord\n";
    }
  }
}else {
  # Parameterized matrix
  #
  # transIFactor = 1 
  #   transVRate = transIRate = mutrate/3
  #
  my $transVRate = $mutRate*($transIFactor/3);
  my $transIRate = ($mutRate - $transVRate)/2;
  my $cpgRate = $cpgFactor*$transVRate;
  $triMatrix = { 
    'AAA' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'AAC' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'AAG' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'AAT' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'ACA' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'ACC' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'ACG' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $cpgRate},
    'ACT' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'AGA' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'AGC' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'AGG' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'AGT' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'ATA' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'ATC' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'ATG' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'ATT' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'CAA' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'CAC' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'CAG' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'CAT' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'CCA' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'CCC' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'CCG' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $cpgRate},
    'CCT' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'CGA' => { 'A' => $cpgRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'CGC' => { 'A' => $cpgRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'CGG' => { 'A' => $cpgRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'CGT' => { 'A' => $cpgRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'CTA' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'CTC' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'CTG' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'CTT' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'GAA' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'GAC' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'GAG' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'GAT' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'GCA' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'GCC' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'GCG' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $cpgRate},
    'GCT' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'GGA' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'GGC' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'GGG' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'GGT' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'GTA' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'GTC' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'GTG' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'GTT' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'TAA' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'TAC' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'TAG' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'TAT' => { 'A' => 0, 'C' => $transVRate, 'G' => $transIRate, 'T' => $transVRate},
    'TCA' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'TCC' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'TCG' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $cpgRate},
    'TCT' => { 'A' => $transVRate, 'C' => 0, 'G' => $transVRate, 'T' => $transIRate},
    'TGA' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'TGC' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'TGG' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'TGT' => { 'A' => $transIRate, 'C' => $transVRate, 'G' => 0, 'T' => $transVRate},
    'TTA' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'TTC' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'TTG' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0},
    'TTT' => { 'A' => $transVRate, 'C' => $transIRate, 'G' => $transVRate, 'T' => 0}
  };


}


# Pre-calculate the sum of the mutation rates for each triplet
my %triMatrixRowSums = ();
foreach my $tri ( keys %{$triMatrix} ) {
  my $rowSum = 0;
  foreach my $base ( 'A', 'C', 'G', 'T' ){
    $rowSum += $triMatrix->{$tri}->{$base};
  }
  $triMatrixRowSums{$tri} = $rowSum;
}
  
# Process remaining options
if ( $options{'verbosity'} ) {
  $verbosity = $options{'verbosity'};
}
my $minFragLen = 100;
if ( $options{'min_frag_len'} ) {
  $minFragLen = $options{'min_frag_len'};
}
my $insertion_rate = 0.08;
if ( $options{'insertion_rate'} ) {
  $insertion_rate = $options{'insertion_rate'};
}
my $deletion_rate = 0.12;
if ( $options{'deletion_rate'} ) {
  $deletion_rate = $options{'deletion_rate'};
}
my $mean_del_size = 1.7;
if ( $options{'mean_del_size'} ) {
  $mean_del_size = $options{'mean_del_size'};
}
my $mean_ins_size = 1.7;
if ( $options{'mean_ins_size'} ) {
  $mean_ins_size = $options{'mean_ins_size'};
}
my $max_del_size = 20;
if ( $options{'max_del_size'} ) {
  $max_del_size = $options{'max_del_size'};
}
my $max_ins_size = 20;
if ( $options{'max_ins_size'} ) {
  $max_ins_size = $options{'max_ins_size'};
}
# Used with 100 seq, 100 generation tree
my $generations_per_unit_time = 5000;
if ( $options{'generations_per_unit_time'} ) {
  $generations_per_unit_time = $options{'generations_per_unit_time'};
}
# Random number generator seed
my $seed;
if ( exists $options{'srand'} ) {
  if ( $options{'srand'} =~ /^\d+$/ ) {
    $seed = $options{'srand'};
    srand($seed);
  }else {
    die "\n\nError: Option -srand must be an integer....received $options{'srand'}\n\n";
  }
}else {
  $seed = srand();
}

#
# Main loop
#
print "TEForwardEvolve.pl - Version $Version\n";
print "======================================\n";
print "Parameters:\n";
print " - random number seed        = $seed\n";
print " - min_fragment_len          = $minFragLen\n";
print " - mutation_rate             = $mutRate\n";
print " - insertion_rate            = $insertion_rate\n";
print " - deletion_rate             = $deletion_rate\n";
print " - mean_ins_size             = $mean_ins_size\n";
print " - max_ins_size              = $max_ins_size\n";
print " - mean_del_size             = $mean_del_size\n";
print " - max_del_size              = $max_del_size\n";
print " - generations_per_unit_time = $generations_per_unit_time\n";
if ( $options{'matrix'} ) {
  print " - custom matrix: $options{'matrix'}\n";
}else {
  if ( $options{'transition_factor'} || $options{'cpg_factor'} ){
    print " - Theoretical matrix: transition_factor = $transIFactor, cpg_factor = $cpgFactor\n";
  }else {
    print " - Theoretical matrix: mutation_CpGx20.txt from the trevover package\n";
  }
}
print "\n\n";

print "Matrix:\n";
foreach my $triplet ( sort { $a cmp $b } keys %{$triMatrix} ) {
  print "$triplet ";
  foreach my $base ( 'A', 'C', 'G', 'T' ) {
    print "\t" . $triMatrix->{$triplet}->{$base};
  }
  print "\n";
}
print "\n\n";


# Perform the simulation
breadthFirstSearch($tref, \&evolveBranch);

# Print multiple alignment
open OUT,">output-msa.fa" or die;
foreach my $node ( @globalNodeList ) {
  if ( $node->{'VALUE'} =~ /int-node/ ) {
    # Do not print internal nodes
  }else {
    print OUT ">" . $node->{'VALUE'} . "\n";
    print OUT "" . $node->{'SEQUENCE'} . "\n";
  }
}
close OUT;

# Print unaligned fasta sequences
open OUT,">output-seqs.fa" or die;
foreach my $node ( @globalNodeList ) {
  if ( $node->{'VALUE'} =~ /int-node/ ) {
    # Do not print internal nodes
  }else {
    print OUT ">" . $node->{'VALUE'} . "\n";
    my $seq = $node->{'SEQUENCE'};
    $seq =~ s/-//g;
    print OUT "$seq\n";
  }
}
close OUT;
 
# TEST
&DFS_avgSubs($tref);
print "Average % Sub Per Extant = " . ($extantSubSum/$extantCnt) . "\n";

# Print final stats
my $subCnt = 0;
my $insCnt = 0;
my $delCnt = 0;
my $sumInsSize = 0;
my $sumDelSize = 0;
foreach my $mut ( @globalCombinedMuts ) {
  if ( $mut->[0] eq "sub" ) {
    $subCnt++;
  }elsif ( $mut->[0] eq "ins" ) {
    $insCnt++;
    $sumInsSize += $mut->[1];
  }elsif ( $mut->[0] eq "del" ) {
    $delCnt++;
    $sumDelSize += $mut->[1];
  }
}
print "Actual Insertion : Deletion Ratio " . ($insCnt/($insCnt+$delCnt)) . " : " . ($delCnt/($insCnt+$delCnt)). "\n";
print "Actual Indel : Substitution Ratio " . ($insCnt/($insCnt+$delCnt+$subCnt)) . " : " . ($subCnt/($insCnt+$delCnt+$subCnt)) . "\n";
print "Actual average insertion length   " . ($sumInsSize/$insCnt) . "\n";
print "Actual average deletion length   " . ($sumDelSize/$delCnt) . "\n";
print "Number of insertion events    $insCnt  ( " . ($insCnt/$subCnt) . " per substitution )\n";
print "Number of deletion events    $delCnt  ( " . ($delCnt/$subCnt ) . " per substitution )\n";
print "Number of substitution events    $subCnt\n";

exit;


######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##

=head2 evolveBranch()

  Use:  my $retVal = evolveBranch( $node, $parent );

    $parameter1:   A generic scalar parameter
    $parameter2:   A generic scalar parameter

  $retVal contains the scalar result of this subroutine.  This
  is a public function and this documentation will print out
  when perldoc is run on this file.

=cut

##-------------------------------------------------------------------------##·
sub evolveBranch{
  my $node = shift;
  my $parent = shift;
  

  # Extract out branch length
  my $branch_length = 0;
  if ( $node->{'VALUE'} =~ /\:([\d\.]+)/ ) {
    $branch_length = $1;
  }

  # If parent is empty this is the root node. This would be a good chance
  # to perform burn-in.  Not sure why that's necessary though.  For now
  # this is unimplemented.
  if ( not defined $parent ) {
    # Do nothing
  }else {
    print "evolveBranch( " . $node->{'VALUE'} . " )\n";
    # Birth new sequence from parent and fragment if necessary
    my $parentSeq = $parent->{'SEQUENCE'};
    my ($parentSeqLen) = ($parentSeq =~ tr/ACGT/ACGT/);
    my $nodeSeq = $parentSeq; 
    # If the node is extant ( ie. doesn't have to be viable ) and
    # we haven't reached our limit and our coin flip is heads
    if ( not exists $node->{'CHILDREN'} ) {
      my $fragLen = int(rand($parentSeqLen-$minFragLen))+$minFragLen;
      my $startPos = int(rand($parentSeqLen-$fragLen));
      my $idx = -1;
      $nodeSeq = "";
      for ( my $i = 0; $i < length($parentSeq); $i++ ) {
        my $b = substr($parentSeq,$i,1);
        $idx++ if ( $b ne "-" );
        if ( $idx >= $startPos && $idx < $startPos + $fragLen ) {
          $nodeSeq .= $b;
        }else {
          $nodeSeq .= "-";
        }
      }
      print " - fragment ( $fragLen, at $startPos ) seq: $nodeSeq\n" if ( $verbosity >= 1 );
    }else {
      print " - full length seq: $nodeSeq\n" if ( $verbosity >= 1 );
    }

    # Calculate total substitution rate
    my $cSeq = $nodeSeq;
    $cSeq =~ s/-//g;
    my $nodeSeqLen = length($cSeq);
    my $total_sub_rate = 0;
    for ( my $i = 0; $i < $nodeSeqLen-2; $i++ ) {
      my $tri = substr($cSeq,$i,3);
      my $rowSum = 0;
      foreach my $subBase ( "A", "C", "G", "T" ) {
        $rowSum += $triMatrix->{$tri}->{$subBase};
      }
      $total_sub_rate += $rowSum;
    }
   
    # Calculate total insertion rate
    # Validate this!
    my $total_insertion_rate = ($insertion_rate*($total_sub_rate/$nodeSeqLen)) * ($nodeSeqLen + 1);

    # Calculate total deletion rate
    my $total_deletion_rate = ($deletion_rate*($total_sub_rate/$nodeSeqLen)) * ($mean_del_size - 1 + $nodeSeqLen);
    
    # Calculate total_mutation_rate 
    my $total_mutation_rate = $total_sub_rate + $total_deletion_rate + $total_insertion_rate;

    # Sample from waiting time distribution with mean = total_mutation_rate and
    # iterate until sum of waiting times exceeds branch length:
    my $waitingTime = -(1/$total_mutation_rate) * log(rand());

    # Calculate generations_in_branch
    my $generations_in_branch = $generations_per_unit_time * $branch_length;

    print " - Branch Length: $branch_length\n" if ( $verbosity >= 2 );
    print " - Generations in branch: $generations_in_branch\n";
    my $elapsedTime = $waitingTime;
    my $subCnt = 0;
    while ( $elapsedTime < $generations_in_branch ) {

      print "  - Waiting time: $waitingTime\n" if ( $verbosity >= 1);
      print "  - Node seq len = $nodeSeqLen\n" if ( $verbosity >= 1);
      print "  - Total mutation rate for node: $total_mutation_rate\n" if ( $verbosity >= 1 );
      if ( $verbosity >= 2 ) {
        print "      * Total substitution rate for node: $total_sub_rate\n";
        print "      * Total insertion rate for node: $total_insertion_rate\n";
        print "      * Total deletion rate for node: $total_deletion_rate\n";
      }

      my $event = rand($total_mutation_rate);
      my $sumRate = 0;
      print "  - Event value = $event\n" if ( $verbosity >= 3);
      #
      # SUBSTITUTION
      #
      if ( ($sumRate += $total_sub_rate) > $event ) {
        my ($newNodeSeq, $new_total_sub_rate, $sub_base_pos, $sub_align_pos, $prev_tri, $new_tri) = 
             &generateSubstitution( $nodeSeq, $total_sub_rate);
        print "  ## Substitution: $prev_tri -> $new_tri @ $sub_align_pos\n" if ( $verbosity >= 1);
        push @globalCombinedMuts, [ "sub", 1 ];
        $nodeSeq = $newNodeSeq;
        $subCnt++;
        $total_sub_rate = $new_total_sub_rate;
      #
      # INSERTION
      #
      }elsif ( ($sumRate += $total_insertion_rate) > $event ) {
        my ($newNodeSeq, $new_total_sub_rate, $ins_base_pos, $ins_align_pos, $insertionSeq, $new_base_len) =
             &generateInsertion( $nodeSeq, $total_sub_rate );
        print "  ## Insertion: @ $ins_align_pos, $insertionSeq\n" if ( $verbosity >= 1);
        $nodeSeq = $newNodeSeq;
        push @globalCombinedMuts, [ "ins", length($insertionSeq) ];
        $total_sub_rate = $new_total_sub_rate;
        $total_insertion_rate = ($insertion_rate*($total_sub_rate/$new_base_len)) * ($new_base_len + 1);
        $total_deletion_rate = ($deletion_rate*($total_sub_rate/$new_base_len)) * ($mean_del_size - 1 + $new_base_len);
      #
      # DELETION
      # 
      }elsif ( ($sumRate += $total_deletion_rate) > $event ) {
        my ($newNodeSeq, $new_total_sub_rate, $del_base_pos, $del_align_pos, $deletionLen, $new_base_len) =
                  &generateDeletion( $nodeSeq, $total_sub_rate, $verbosity );
        print "  ## Deletion: @ $del_align_pos, length = $deletionLen\n" if ( $verbosity >= 1);
        $nodeSeq = $newNodeSeq;
        push @globalCombinedMuts, [ "del", $deletionLen ];
        $total_sub_rate = $new_total_sub_rate;
        $total_insertion_rate = ($insertion_rate*($total_sub_rate/$new_base_len)) * ($new_base_len + 1);
        $total_deletion_rate = ($deletion_rate*($total_sub_rate/$new_base_len)) * ($mean_del_size - 1 + $new_base_len);
      }
      print "NEW SEQ = $nodeSeq\n" if ( $verbosity >= 4 );

      # re-Calculate total_mutation_rate 
      $total_mutation_rate = $total_sub_rate + $total_deletion_rate + $total_insertion_rate;
      print "  - Adjusted total mutation rate for node: $total_mutation_rate\n" if ( $verbosity >= 1 );

      $waitingTime = -(1/$total_mutation_rate) * log(rand());
      $elapsedTime += $waitingTime;
    }
    $node->{'SEQUENCE'} = $nodeSeq;
    $node->{'SUBCNT'} = $subCnt;
    $node->{'CUMSUBCNT'} = $subCnt;
    if ( exists $parent->{'CUMSUBCNT'} && defined $parent->{'CUMSUBCNT'} ) {
      $node->{'CUMSUBCNT'} += $parent->{'CUMSUBCNT'};
    }
  } # If root node
}


##-------------------------------------------------------------------------##

=head2 chooseIndelLength()

  Use:  my $retVal = chooseIndelLength( $exponent, $max );

    $exponent:    The exponent of the distribution
    $max     :    The max value to accept from the distribution 

  $retVal the length of an indel chosen from a power law probability
  distribution.

  Adapted from INDELIBLE, which in turn pulled this algorithm from 
  DAWG (Cartwright, 2005).  
 
  They draw from a Zipfian distribution with a parameter 
  $exponent > 1.0. Ref: Devroye Luc (1986) Non-uniform random 
  variate generation. Springer-Verlag: Berlin. p551
 
  If this ends up being too slow INDELIBLE used the rejection-inversion 
  method of Hormann and Derflinger (1996) to perform the same sample.

=cut

##-------------------------------------------------------------------------##·
sub chooseIndelLength{
  my $exponent = shift;
  my $max = shift;

  my $x;
  do { 
    my $b = 2.0**($exponent-1.0);
    $x = 0;
    my $t = 0;
    do { 
      $x = floor(rand()**(-1.0/$exponent-1.0));
      $t = (1.0+1.0/$x)**($exponent-1.0);
    }while( rand()*$x*($t-1.0)*$b > $t*($b-1.0));
  }while ( int($x) > $max );
  return int($x);
}


#-------------------------------------------------------------------------##

=head2 generateDeletion()

  Use:  my ($newNodeSeq, $new_total_sub_rate, $del_base_pos,
            $del_align_pos, $deletionLen, $new_base_len ) = 
                      generateDeletion( $nodeSeq, $total_sub_rate, 
                                        $verbosity );

    $nodeSeq         :    The exponent of the distribution
    $total_sub_rate  :    The max value to accept from the distribution 
    $verbosity       :

  Returns.
  Uses power-law 
  Currently warns if sequence length is less than 4bp at the start.
  Returns without deletion if sequence is zero bp.
  Deletions that run off the end of the sequence ( because start
  position is picked independent of size ) are equivalent to a deletion
  of size (sequence_length - start_position).
  

    $newNodeSeq         :
    $new_total_sub_rate :
    $del_base_pos       :
    $del_align_pos      :
    $deletionLen        :
    $new_base_len       :

=cut

##-------------------------------------------------------------------------##·
sub generateDeletion{
  my $nodeSeq = shift;
  my $total_sub_rate = shift;
  my $verbosity = shift;

  # Create a gapless sequence from the aligned sequence and
  # generate a set of indices back to the gapped sequence.
  # 
  #  E.g             aligned sequence: A--ACCA---TA---
  #           gapless sequence (cSeq): A  A  C  C  A  T  A
  #         source indices absPosList: 0  3  4  5  6 10 11
  my @absPosList = ();
  my $cSeq = "";
  for ( my $i = 0; $i < length($nodeSeq); $i++ ) {
    my $base = substr($nodeSeq, $i, 1);
    if ( $base ne "-" ) {
      $cSeq .= $base;
      push @absPosList, $i;
    }
  }

  # Perhaps too verbose...but could help flag a parameter issue.
  if ( length($cSeq) < 4 ) {
    warn "WARNING: generateDeletion - sequence length has dropped below 4bp!\n";
  }

  # Nothing left of this sequence!  Hopefully we won't get called in
  # such a case...but if so
  if ( length($cSeq) == 0 ) {
    return ($nodeSeq, $total_sub_rate, 0, 0, 0, 0);
  }

  # Pick location to start deleting (zero based)
  my $del_base_pos = int(rand(length($cSeq)));
  my $del_align_pos = $absPosList[$del_base_pos];
  if ( $del_align_pos eq "" ) {
     die "ERROR: generateDeletion - del_base_pos of $del_base_pos doesn't have an alignment position!\n";
  }

  # Currently uses a power-law distribution for picking deletion sizes
  my $deletionLen = chooseIndelLength($mean_del_size, $max_del_size);
  #   NOTE: Location/Size of the deletion is not restricted by sequence edges.
  #         If a deletion will run off the sequence the mutation is still 
  #         performed abeit with a smaller effect. It is recorded as the 
  #         effective length for accounting.
  $deletionLen = length($cSeq) - $del_base_pos if ( $del_base_pos + $deletionLen > length($cSeq) ); 
  my $del_align_endpos = $absPosList[$del_base_pos + $deletionLen - 1];
  if ( $del_align_endpos eq "" ) {
     die "ERROR: generateDeletion - del_base_endpos of " . ($del_base_pos + $deletionLen - 1) . " doesn't have an alignment position!\n";
  }
  if ( $verbosity >= 3 ) {
    print "INFO: generateDeletion - deleting $deletionLen" . "bp starting at alignment_pos=$del_align_pos/$del_base_pos\n";
  }
  
  # Adjust total substitution rate due to composition changes
  # - Obtain at most, two left flanking bases.  Fewer if we
  #   are close to the start of the sequnce.
  my $new_total_sub_rate = $total_sub_rate;
  my $left_seq = "";
  if ( $del_base_pos == 1 ) {
    # Only one flanking base
    $left_seq = substr($cSeq,$del_base_pos-1,1);
  }elsif ( $del_base_pos > 1 ) {
    # Two flanking bases
    $left_seq = substr($cSeq,$del_base_pos-2,2);
  } 

  # - Obtain at most, two right flanking bases.  Fewer if we
  #   are close to the end of the sequnce.
  my $right_seq = "";
  my $del_base_endpos = ( $del_base_pos + $deletionLen - 1 );
  if ( $del_base_endpos == length($cSeq)-1 ) {
    # One flanking base
    $right_seq = substr($cSeq,$del_base_endpos+1,1);
  }elsif ( $del_base_endpos < length($cSeq)-1 ) {
    # Two flanking bases
    $right_seq = substr($cSeq,$del_base_endpos+1,2);
  } 

  print "INFO: generateDeletion - deletion context <$left_seq>" .  substr($cSeq,$del_base_pos,$deletionLen) . "<$right_seq>\n"
    if ( $verbosity >= 4 );

  my $combined = $left_seq . substr($cSeq,$del_base_pos,$deletionLen) . $right_seq;
  for ( my $i = 0; $i < length($combined)-2; $i++ )
  {
    print "INFO: generateDeletion - removing triplet rate for " . substr($combined,$i,3) . "\n"
      if ( $verbosity >= 4 );
    $new_total_sub_rate -= $triMatrixRowSums{substr($combined,$i,3)};
  }
  $combined = $left_seq . $right_seq;
  for ( my $i = 0; $i < length($combined)-2; $i++ )
  {
    print "INFO: generateDeletion - adding triplet rate for " . substr($combined,$i,3) . "\n"
      if ( $verbosity >= 4 );
    $new_total_sub_rate += $triMatrixRowSums{substr($combined,$i,3)};
  }

  # Perform operation on aligned sequence
  my $newNodeSeq = $nodeSeq;
  substr($newNodeSeq,$del_align_pos,$del_align_endpos - $del_align_pos + 1) = '-'x($del_align_endpos - $del_align_pos + 1);

  my $new_base_len = length($cSeq) - $deletionLen;

  die "ERROR: generateDeletion - aligned sequence length changed during deletion!\n" if ( length($newNodeSeq) != length($nodeSeq) );
  
  # NOTE: Caller must update total_insertion_rate and total_deletion_rate for 
  #       changes in sequence length
  return ($newNodeSeq, $new_total_sub_rate, $del_base_pos, $del_align_pos, $deletionLen, $new_base_len);
}


#-------------------------------------------------------------------------##

=head2 generateInsertion()

  Use:  my ($newNodeSeq, $new_total_sub_rate, $ins_base_pos,
            $ins_align_pos, $insertionSeq, $new_base_len ) = 
                      generateInsertion( $nodeSeq, $total_sub_rate, 
                                        $verbosity );

    $nodeSeq         :    The current aligned sequence for the node
    $total_sub_rate  :    The current total substitution rate for the node
    $verbosity       :    The level of verbosity for logging

  Returns.
  Uses power-law 
  Currently warns if sequence length is less than 4bp at the start.
  Returns without insertion if sequence is zero bp.
  Insertions can occur before/after the sequence.
  
    $newNodeSeq         :
    $new_total_sub_rate :
    $ins_base_pos       :
    $ins_align_pos      :
    $insertionSeq       :
    $new_base_len       :

=cut

##-------------------------------------------------------------------------##·
sub generateInsertion{
  my $nodeSeq = shift;
  my $total_sub_rate = shift;
  my $verbosity = shift;

  # Create a gapless sequence from the aligned sequence and
  # generate a set of indices back to the gapped sequence.
  # 
  #  E.g             aligned sequence: A--ACCA---TA---
  #           gapless sequence (cSeq): A  A  C  C  A  T  A
  #         source indices absPosList: 0  3  4  5  6 10 11
  my @absPosList = ();
  my $cSeq = "";
  for ( my $i = 0; $i < length($nodeSeq); $i++ ) {
    my $base = substr($nodeSeq, $i, 1);
    if ( $base ne "-" ) {
      $cSeq .= $base;
      push @absPosList, $i;
    }
  }
  # Insertions occur prior to the start position selected.  In order 
  # to support insertions before *and* after the sequence it is necessary
  # to reference a sequence position *after* the last base.  This next
  # line simply creates this out-of-sequence lookup.  For example:
  #
  #    sequence   =  AAACA----
  #    cSeq       =  AAACA
  #    absPosList =  012349
  #
  # If we decide to insert "TA" the end of this sequence the new sequence
  # would show up right shifted after the last base as:
  #
  #   new sequence = AAACA----TA
  #
  push @absPosList, length($nodeSeq);

  # Perhaps too verbose...but could help flag a parameter issue.
  if ( length($cSeq) < 4 ) {
    warn "WARNING: generateInsertion - sequence length has dropped below 4bp!\n";
  }
   
  # Nothing left of this sequence!  There is no real reference for
  # where to insert in this case nor is it desirable to simulate the
  # birth of something new here.  Hopefully we won't get called in
  # such a case...but if so
  if ( length($cSeq) == 0 ) {
    return ($nodeSeq, $total_sub_rate, 0, 0, 0, 0);
  }

  # Pick location ( insertions occur to the left of the position selected )
  my $ins_base_pos = rand(length($cSeq)+1);
  my $ins_align_pos = $absPosList[$ins_base_pos];

  my $insertionLen = chooseIndelLength($mean_ins_size, $max_ins_size);
  my $insertionSeq = "";
  ## TODO: Should have a parameter for the background base frequencies for insertions
  my @bases = ('A','C','G','T');
  for ( my $i = 0; $i < $insertionLen; $i++ ){
    $insertionSeq .= $bases[int(rand(4))];
  }
  #print "INFO: generateInsertion - inserting $insertionSeq ($insertionLen" . "bp) starting prior to alignment_pos=$ins_align_pos/$ins_base_pos\n";
  
  # Adjust total substitution rate due to composition changes
  # - Obtain at most, two left flanking bases.  Fewer if we
  #   are close to the start of the sequnce.
  my $new_total_sub_rate = $total_sub_rate;
  my $left_seq = "";
  if ( $ins_base_pos == 1 ) {
    $left_seq = substr($cSeq,$ins_base_pos-1,1);
  }elsif ( $ins_base_pos > 1 ) {
    $left_seq = substr($cSeq,$ins_base_pos-2,2);
  } 
  # - Obtain at most, two right flanking bases.  Fewer if we
  #   are close to the end of the sequnce.
  my $right_seq = "";
  if ( $ins_base_pos == length($cSeq)-1 ) {
    $right_seq = substr($cSeq,$ins_base_pos,1);
  }elsif ( $ins_base_pos < length($cSeq)-1 ) {
    $right_seq = substr($cSeq,$ins_base_pos,2);
  } 
  #print "INFO: generateInsertion - insertion context <$left_seq><$right_seq>\n";
  my $combined = $left_seq . $right_seq;
  for ( my $i = 0; $i < length($combined)-2; $i++ )
  {
    #print "INFO: generateInsertion - removing triplet rate for " . substr($combined,$i,3) . "\n";
    $new_total_sub_rate -= $triMatrixRowSums{substr($combined,$i,3)};
  }
  $combined = $left_seq . $insertionSeq . $right_seq;
  for ( my $i = 0; $i < length($combined)-2; $i++ )
  {
    #print "INFO: generateInsertion - adding triplet rate for " . substr($combined,$i,3) . "\n";
    $new_total_sub_rate += $triMatrixRowSums{substr($combined,$i,3)};
  }

  # Perform operation on aligned sequence
  my $newNodeSeq = $nodeSeq;
  substr($newNodeSeq,$ins_align_pos,0) = $insertionSeq;

  die "ERROR: generateInsertion - aligned sequence length changed by incorrect amount during insertion!\n" 
       if ( length($newNodeSeq) != length($nodeSeq)+$insertionLen );

  #
  # Pad all previous sequences to accomodate new insertion
  #
  my $newAlignedLen = length($newNodeSeq);
  foreach my $node ( @globalNodeList ) {
    if ( exists $node->{'SEQUENCE'} ) {
      substr($node->{'SEQUENCE'}, $ins_align_pos,0) = '-'x($insertionLen);
      #die "ERROR: generateInsertion - aligned sequence length changed by incorrect amount during insertion!\n" 
      #   if ( length($node->{'SEQUENCE'}) != $newAlignedLen );
    }
  }

  my $new_base_len = length($cSeq) + $insertionLen;

  # NOTE: Caller must update total_insertion_rate and total_deletion_rate for 
  #       changes in base length
  return ($newNodeSeq, $new_total_sub_rate, $ins_base_pos, $ins_align_pos, $insertionSeq, $new_base_len);
  
}


#-------------------------------------------------------------------------##

=head2 generateSubstitution()

  Use:  my ($newNodeSeq, $new_total_sub_rate, $sub_base_pos,
            $sub_align_pos, $prev_tri, $new_tri ) = 
                      generateInsertion( $nodeSeq, $total_sub_rate, 
                                        $verbosity );

    $nodeSeq         :    The current aligned sequence for the node
    $total_sub_rate  :    The current total substitution rate for the node
    $verbosity       :    The level of verbosity for logging

  Returns.
  
    $newNodeSeq         :
    $new_total_sub_rate :
    $sub_base_pos       :
    $sub_align_pos      :
    $prev_tri           :
    $new_tri            :

=cut

##-------------------------------------------------------------------------##·
sub generateSubstitution{
  my $nodeSeq = shift;
  my $total_sub_rate = shift;
  my $verbosity = shift;

  # Create a gapless sequence from the aligned sequence and
  # generate a set of indices back to the gapped sequence.
  # 
  #  E.g             aligned sequence: A--ACCA---TA---
  #           gapless sequence (cSeq): A  A  C  C  A  T  A
  #         source indices absPosList: 0  3  4  5  6 10 11
  my @absPosList = ();
  my $cSeq = "";
  for ( my $i = 0; $i < length($nodeSeq); $i++ ) {
    my $base = substr($nodeSeq, $i, 1);
    if ( $base ne "-" ) {
      $cSeq .= $base;
      push @absPosList, $i;
    }
  }
  
  #
  # Based on an approach used in trevolver:
  #
  # Using the sum of mutation rates (M) for all sites within
  # the sequence, randomly choose a value between 0 and M to
  # use as a threshold for identifying the location to 
  # mutate.  E.g
  #
  #   rateThresh = rand(total_sub_rate) = 2.33E-05
  #
  #   ACCGTCGATAG
  #   ACC          [+7.5E-07]
  #    CCG         [+5.0E-06] = 6.26E-06
  #     CGT        [+5.5E-06] = 1.18E-05
  #      GTC       [+7.5E-07] = 1.25E-05
  #       TCG      [+5.5E-06] = 1.80E-05
  #        CGA     [+5.5E-06] = 2.35E-05 This site exceeds
  #                                      2.33E-05
  #
  #     CGA->CAA [+5.0E-06] = 2.30E-05
  #     CGA->CCA [+2.5E-07] = 2.33E-05
  #     CGA->CGA [+0]       = 2.23E-05
  #     CGA->CTA [+2.5E-07] = 2.35E-05 T mutation selected.
  #
  #print "sub rate thresh: $rateThresh / $total_sub_rate\n";
  my $sum_sub_rates = 0;
  my $sub_base_pos = -1;
  my $sub_align_pos = -1;
  my $prev_tri = "";
  my $new_tri = "";
  my $left_prev_tri = "";
  my $left_new_tri = "";
  my $right_prev_tri = "";
  my $right_new_tri = "";
  my $new_total_sub_rate = $total_sub_rate;
  my $sub_base = "";
  my $rateThresh;
  # NOTES: 
  #   - First and last base cannot be mutated with this method.
  #   - Floating point precision errors compound in the calculation 
  #     of total_sub_rate. This can lead to rare cases where the
  #     subsitution is not located within the sequence.  In these
  #     cases we draw another random sample.
  while ( $sub_base eq "" ) {
    $rateThresh = rand($total_sub_rate);
    for ( my $i = 0; $i < length($cSeq)-2; $i++ ) {
      my $tri = substr($cSeq,$i,3);
      if ( $sum_sub_rates + $triMatrixRowSums{$tri} > $rateThresh  ) {
        # Found position...find mutation
        foreach my $base ( 'A', 'C', 'G', 'T' ){
          $sum_sub_rates += $triMatrix->{$tri}->{$base};
          if ( $sum_sub_rates > $rateThresh ) {
            $sub_base = $base;
            $sub_base_pos = $i+1;
            $sub_align_pos = $absPosList[$i+1];
            #print "Found mutation to perform: $tri -> $sub_base, base_pos=$sub_base_pos aligned_pos=$sub_align_pos\n";
            #
            # Center tripplet e.g. flank mutation flank
            #
            $prev_tri = $tri;
            $new_tri = $tri;
            substr($new_tri,1,1) = $sub_base;
            # Adjust total substitution rate due to composition changes
            $new_total_sub_rate -= $triMatrixRowSums{$prev_tri};
            $new_total_sub_rate += $triMatrixRowSums{$new_tri};
    
            #
            # Left triplet e.g. flank flank mutation
            #
            if ( $i > 0 ) {
              $left_prev_tri = substr($cSeq,$i-1,3);
              $left_new_tri = $left_prev_tri;
              substr($left_new_tri,2,1) = $sub_base;
              # Adjust total substitution rate due to composition changes
              $new_total_sub_rate -= $triMatrixRowSums{$left_prev_tri};
              $new_total_sub_rate += $triMatrixRowSums{$left_new_tri};
            }
            
            #
            # Right triplet e.g. mutation flank flank
            # 
            if ( $i < length($cSeq)-3 ) {
              $right_prev_tri = substr($cSeq,$i+1,3);
              $right_new_tri = $right_prev_tri;
              substr($right_new_tri,0,1) = $sub_base;
              # Adjust total substitution rate due to composition changes
              $new_total_sub_rate -= $triMatrixRowSums{$right_prev_tri};
              $new_total_sub_rate += $triMatrixRowSums{$right_new_tri};
            }
            last;
          }
        } # foreach base
      }else {
        # site not higher than threshold
        $sum_sub_rates += $triMatrixRowSums{$tri};
      }
      # Site found!
      last if ( $sub_base );
    }
    if ( ! $sub_base && $verbosity >= 2 ) {
      print "  !! Lack of floating point precision caused rate_threshold to be set too high.\n";
    }
  }

  my $newNodeSeq = $nodeSeq;
  substr($newNodeSeq,$sub_align_pos,1) = $sub_base;
  if ( substr($nodeSeq,$sub_align_pos,1) eq "-" ||  length($newNodeSeq) != length($nodeSeq) )
  {
    die "Substitution ERROR: origAlignLen=".length($nodeSeq)." subAlignPos=$sub_align_pos origBase=" . substr($nodeSeq,$sub_align_pos,1) ."\n";
  }

  return ($newNodeSeq, $new_total_sub_rate, $sub_base_pos, $sub_align_pos, $prev_tri, $new_tri);
}


#-------------------------------------------------------------------------##

=head2 fixMissingNodeNames()

  Use: &fixMissingNodeNames($node, $parent);

    $node     :    ...
    $parent   :    

  Modifies the tree datastructure...

=cut

##-------------------------------------------------------------------------##·
sub fixMissingNodeNames {
  my $node = shift;
  my $parent = shift;

  my $idx = 1;
  my $nodeName = $node->{'VALUE'};
  # s1:0.1 or :0.1
  if ( $nodeName =~ /([^;:]*)\:([\d\.]+)/ ) {
    my $distance = $2;
    if ( $1 eq "" ) {
      while ( exists $globalNodeNamesHash{"Unnamed".$idx} ) {
        $idx++;
      }
      $globalNodeNamesHash{"Unnamed".$idx}++;
      #print "Fixing nodename >$nodeName< to Unnamed$idx\n";
      $node->{'VALUE'} = "Unnamed$idx:$distance";
    }
  # s1 or ;
  }elsif ( $nodeName =~ /^([^:]+)$/ ){
    if ( $1 eq ";" || $1 eq "" ) {
      while ( exists $globalNodeNamesHash{"Unnamed".$idx} ) {
        $idx++;
      }
      $globalNodeNamesHash{"Unnamed".$idx}++;
      #print "Fixing nodename2 >$nodeName< to Unnamed$idx\n";
      $node->{'VALUE'} = "Unnamed$idx";
    }
  }else {
    print "Didn't expect this nodename: $nodeName\n";
  }
  push @globalNodeList, $node;
}


## TODO: This is a testing routine...and based on Tigger1 length
sub DFS_avgSubs {
  my $f = shift;

  if ( exists $f->{'CHILDREN'} && @{$f->{'CHILDREN'}} ) {
    for my $node ( @{$f->{'CHILDREN'}} ) {
      &DFS_avgSubs($node);
    }
  }else {
    $extantSubSum += $f->{'CUMSUBCNT'}/2418;
    $extantCnt++;
  }
}


##-------------------------------------------------------------------------##

=head2 breadthFirstSearch()

  Use: &breadthFirstSearch($root, $operation);

    $root        :    ...
    $operation   :    

  Basic recursive breadth first search algorithm

=cut

##-------------------------------------------------------------------------##·
sub breadthFirstSearch{
  my $root = shift;
  my $operation = shift;

  # Start of search, treat 'seen' as flip flop
  my $seenState = 1;
  if ( exists $root->{'seen'} ) {
    if ( $root->{'seen'} == 1 ) {
      $seenState = 0;
    }
  }

  my @queue = ();
  
  # Start with the root node
  push @queue, $root;

  # Do something with the root node here:
  $operation->( $root, undef );
  $root->{'seen'} = $seenState;

  # Continue with the children
  while ( @queue ) {
    my $node = shift @queue;
    my @children = ();
    if ( exists $node->{'CHILDREN'} ) {
      foreach my $child ( @{$node->{'CHILDREN'}} ) {
        push @children, $child;
      }
    }
    foreach my $child ( @children ) {
      if ( ! exists $child->{'seen'} || $child->{'seen'} != $seenState ){
        $child->{'seen'} = $seenState;
        # Begin : Do something with the child node here:
        $operation->( $child, $node );
        # End :  Do something with child node 
        push @queue, $child;
      }
    }
  }
}

1;
