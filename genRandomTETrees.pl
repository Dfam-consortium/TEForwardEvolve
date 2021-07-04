#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) program_name
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      This is a template for generic perl scripts.  It
##      includes an area for POD documentation (ie. perldoc this-file )
##      and a generic way to store startup parameters in a file
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

program_name - some description

=head1 SYNOPSIS

  program_name [-version]

=head1 DESCRIPTION

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

my $Version = "0.1";

my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-version', # print out the version and exit
    '-dna',
    '-line'
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

unless ( $options{'dna'} || $options{'line'} ) {
  die "Must specify either -dna or -line\n";
}

my $nw;

if ( $options{'dna'} ) {
  #
  # DNA Transposon
  #   - Many generations of non-competent copies from a single live sequence, bush-like phylogeny
  #   - The parameter "totalGenerations" controls the amount of time post-extinction represented
  #     in the extant sequence branch lengths.
  #
  my $extantNodeGoal = 100;
  my $avgGenerationsBetweenParents = 10;
  my $totalGenerations =  100;
  
  my $numExtantNodes = 1;
  my $root = { 'VALUE' => "node-0",
               'TOROOT' => 0,
               'PARENT' => undef,
               'CHILDREN' => [] };
  
  my @nodes;
  push @nodes, $root;
  
  # Simulate the active birth of sequences
  my $idx = 1;
  while ( $numExtantNodes < $extantNodeGoal ) {
    my $parent;
    if ( $#nodes == 0 ) {
      $parent = $nodes[0];
    }else {
      $parent = $nodes[int(rand($#nodes))];
    }
    #print "nodes = " . scalar(@nodes) . " and parent choice is $parent\n";
    my $branchLen = int(rand($avgGenerationsBetweenParents+1));
  
    if ( @{$parent->{'CHILDREN'}} > 0 ) {
      $numExtantNodes++;
    }
  
    my $child = { 'VALUE' => "node-$idx",
               'TOROOT' => $branchLen + $parent->{'TOROOT'},
               'PARENT' => $parent,
               'CHILDREN' => [] };
    $idx++;
    push @{$parent->{'CHILDREN'}}, $child;
    push @nodes, $child;
  }
  
  # Adjust extant branch lengths to reach the target totalGenerations.
  my $changes;
  do { 
    $changes = 0;
    for ( my $i = 1; $i <= $#nodes; $i++ ) {
      if ( @{$nodes[$i]->{'CHILDREN'}} == 0 && @{$nodes[$i]->{'PARENT'}->{'CHILDREN'}} == 1 ) {
        $nodes[$i]->{'PARENT'}->{'CHILDREN'} = [];
        $changes++;
      }
    }
  } while ( $changes > 0);
   
  for ( my $i = 0; $i <= $#nodes; $i++ ) {
    if ( @{$nodes[$i]->{'CHILDREN'}} == 0 ) {
      $nodes[$i]->{'VALUE'} .= ":" . ( $totalGenerations - $nodes[$i]->{'PARENT'}->{'TOROOT'} );
    }else {
      if ( not defined $nodes[$i]->{'PARENT'}->{'TOROOT'} ) {
        $nodes[$i]->{'VALUE'} .= ":" . ($nodes[$i]->{'TOROOT'});
      }else { 
        $nodes[$i]->{'VALUE'} .= ":" . ($nodes[$i]->{'TOROOT'} - $nodes[$i]->{'PARENT'}->{'TOROOT'});
      }
    }
  }
  
  #        A            ((E, F)B, C, D)A
  #      / | \
  #     B  C  D      - 
  #    /\              
  #   E  F
  #
  #print "Tree: " . Dumper($root) . "\n";
  $nw = "";
  &dfs($root);
  $nw =~ s/\(,/\(/g;
  $nw =~ s/^,//g;

  print "DNA Transposon Random Tree:\n";
  print "  - extantNodeGoal = $extantNodeGoal\n";
  print "  - avgGenerationsBetweenParents = $avgGenerationsBetweenParents\n";
  print "  - totalGenerations = $totalGenerations\n";
  print "Newick: $nw\n";

}else {
  #
  # LINE-Like Phylogeny
  #
  my $totalExtantNodes = 100;
  # Produce a distribution of children between 2-#
  my $avgChildrenPerSeq = 5; # must be > 2 
  my $avgGenerationsBetweenParents = 4.7;
  my $totalGenerations =  100;
  
  
  my $numNodes = 0;
  my $root = { 'VALUE' => "node-0",
               'CHILDREN' => [] };
  my $parent = $root;
  my $sumGenerations = 0;
  my $numParentGens = 0;
  my $lastFlag = 0;
  do { 
       # How many children?
       my $numChildren = int(rand($avgChildrenPerSeq-1)+2);
       if ( $totalExtantNodes - $numNodes < 4 ) {
         $numChildren = 3;
         $lastFlag = 1;
       }
       print "num children = $numChildren\n";
       if ( $numChildren ) {
         my $nextParent;
         for ( my $i = 0; $i < $numChildren; $i++ ) {
           my $node;
           if ( $i == 0 && $totalExtantNodes - ( $numNodes + $numChildren) >= 2 ) 
           {
             # New parent
             $numParentGens = int(rand($avgGenerationsBetweenParents+1));
             $node = { 'VALUE' => "node-$numNodes:" . ( $numParentGens ) . ".0",
                       'CHILDREN' => [] };
             $nextParent = $node;
           }else { 
             # Extant child....scale to total generation length
             $node = { 'VALUE' => "node-$numNodes:" . ( $totalGenerations - $sumGenerations ) . ".0",
                       'CHILDREN' => [] };
             $numNodes++;
           }
           push @{$parent->{'CHILDREN'}}, $node;
         }
         $parent = $nextParent;
         $sumGenerations += $numParentGens;
       }
  }while ( $numNodes < $totalExtantNodes );
  
  
  #        A            ((E, F)B, C, D)A
  #      / | \
  #     B  C  D      - 
  #    /\              
  #   E  F
  #
  #print "Tree: " . Dumper($root) . "\n";
  $nw = "";
  &dfs($root);
  $nw =~ s/\(,/\(/g;
  $nw =~ s/^,//g;

  print "LINE-like Random Tree:\n";
  print "   - totalExtantNodes = $totalExtantNodes\n";
  print "   - avgChildrenPerSeq = $avgChildrenPerSeq\n";
  print "   - avgGenerationsBetweenParents = $avgGenerationsBetweenParents\n";
  print "   - totalGenerations = $totalGenerations\n";
  print "Newick: $nw\n\n";

  my $verify = $nw;
  $verify =~ s/int-node//g;
  my ($cnt) = ($verify =~ s/node/node/g);
  print "Total extants = $cnt\n";
  print "numNodes = $numNodes\n";
}


######################## S U B R O U T I N E S ############################


##-------------------------------------------------------------------------##

=head2 dfs()

  Use:  &dfs( $tree ); 

      $tree = { 'VALUE' => 'node1:0.1',
                'CHILDREN'} => [ { 'VALUE' => 'node2:0.2',
                                   'CHILDREN' => [] } ]
              }

  Perform a depth-first-search of the tree ( recursively ) and
  write Newick(ish) like data to the global $nw variable.  Afterwards
  two transformations need to be performed to have "true" Newick:

  $nw =~ s/\(,/\(/g;
  $nw =~ s/^,//g;

=cut

##-------------------------------------------------------------------------##·
sub dfs {
  my $f = shift;
  
  if ( @{$f->{'CHILDREN'}} ) {
    $nw .= ",(";
    for my $node ( @{$f->{'CHILDREN'}} ) {
      &dfs($node);
    }
    $nw .= ")" . "int-" . $f->{"VALUE"};
  }else {
    $nw .= "," . $f->{"VALUE"};
  }
}

1;
