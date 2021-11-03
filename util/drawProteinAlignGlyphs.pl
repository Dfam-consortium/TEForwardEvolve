#!/usr/local/bin/perl
use strict;
use warnings;
use SVG;
use Data::Dumper;

my $showUnion = 1;
#my $absScale = 1;
my @methodOrder = ("refiner", "clustalo", "muscle", "probcons", "mafft", "dialign", "kalign");

print "<html>\n";
my %colors = ( "clustalo" => "rgb(161, 34, 148)",
               "dialign" => "rgb(237, 19, 7)",
               "refiner" => "rgb(217, 214, 52)",
               "muscle" => "rgb(81, 0, 255)",
               "mafft" => "rgb(33, 145, 50)",
               "kalign" => "rgb(235, 162, 215)",
               "probcons" => "rgb(255, 128, 0)" );

my %useProts = ( "Arthur2" => "KAF0768179.1",
                 "Tigger10" => "XP_005302607.1",
                 "Zaphod" => "KAE9522890.1",
                 "Zaphod2" => "XP_022019873.1",
               );
my %famData = ();
while ( <> ){
  # Zaphod2	846	refiner	KAE9522890.1	5e-06	[451,489]
  my @flds = split();
  my $protein = $flds[3];
  my $method = $flds[2];
  my $protLen = $flds[1];
  my $family = $flds[0];
  my $evalue = $flds[4];

  # FILTER
  if ( $protein =~ /\#/ || $useProts{$family} eq $protein ) {
  }else { next; }
  next if ( $evalue > 0.001 );

  $famData{$family} = {} if ( ! exists $famData{$family} );
  $famData{$family}->{$protein} = {} if ( ! exists $famData{$family}->{$protein} );
  
  my $record = $famData{$family}->{$protein};
  $record->{'length'} = $protLen;
  $record->{'annotData'} = {} if ( ! exists $record->{'annotData'} );
  $record->{'annotData'}->{$method} = [] if ( ! exists $record->{'annotData'}->{$method} );
  my @an = ();
  for ( my $i = 5; $i <= $#flds; $i++ ) {
    if ( $flds[$i] =~ /\[\s*(\d+)\s*,\s*(\d+)\s*\]/ ) {
      push @an, [$1,$2];
    }
  }
  if ( $showUnion ) {
    if ( @{$record->{'annotData'}->{$method}} == 1 ) {
      push @{$record->{'annotData'}->{$method}->[0]}, @an;
    }else {
      push @{$record->{'annotData'}->{$method}}, \@an;
    }
  }else {
    push @{$record->{'annotData'}->{$method}}, \@an;
  }
}
#print "Dumper: " . Dumper(\%famData) . "\n";
#exit;

foreach my $family ( keys(%famData) ) {
  print "<h2>$family</h2>\n";
print "<div>\n";
foreach my $method ( keys %colors ) {
  my $color = $colors{$method};
  print "<div>\n";
  print "<svg width='10' height='10' xmlns='http://www.w3.org/2000/svg' xmlns:svg='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n";
  print "<rect width='10' height='10' style='fill:" . $color . 
        ";stroke:" . $color . "'/></svg>\n";
  print "$method</div>\n";
}
print "</div>\n";


  foreach my $prot ( keys(%{$famData{$family}} ) ) {
    my $pRec = $famData{$family}->{$prot};
    print "<h3>$prot</h3>\n";
    print "<i><b>Protein Length = $pRec->{'length'}</b></i><br>\n";
    drawSVG($pRec->{'length'}, $pRec->{'annotData'}, $prot, $family);
  }
  print "<hr>\n";
}
  
  
  
 
    
#my $proteinLen = 846;
#my $annotData = { "clustalo" =>  [ [ [823,844] ] ],
#                  "dialign" => [ [ [803,829,1], [831,846,1] ] ],
#                  "refiner" => [ [ [801,846,2], [678,742,2], [756,801,2] ],
#                                 [ [506,576,2], [439,505,2] ],
#                                 [ [423,498,2] ] ],
#                  "mafft" =>    [ [ [383,427,3] ] ],
#                  "probcons" => [ [ [450,472,4] ] ]
#                 };
##Zaphod	clustalo	KAE9522890.1	4e-05	[823,844]
#Zaphod	dialign	KAE9522890.1	6e-06	[803,829]	[831,846]
#Zaphod	refiner	KAE9522890.1	7e-39	[801,846]	[678,742]	[756,801]
#Zaphod	refiner	KAE9522890.1	5e-34	[506,576]	[439,505]
#Zaphod	refiner	KAE9522890.1	4e-08	[423,498]
#Zaphod	mafft	KAE9522890.1	5e-05	[383,427]
#Zaphod	probcons	KAE9522890.1	4e-04	[450,472]



 

sub drawSVG {
  my $proteinLen = shift;
  my $annotData = shift;
  my $protName = shift;
  my $famName = shift;

# create an SVG object
my $svg= SVG->new( width => 200, height => 200);
 
my $gradient = $svg->gradient(
    -type => "linear",
    id    => "gradient_1",
    x1=>"0%", y1=>"0%", x2=>"0%", y2=>"100%"
);

$gradient->stop(
    id     => "stop_1",
    offset  => "0%",
    'stop-color'  => "rgb(85, 145, 201)",
);
$gradient->stop(
    id     => "stop_2",
    offset  => "50%",
    'stop-color'  => "rgb(85, 145, 201)",
    'stop-opacity' => "0"
);
$gradient->stop(
    id     => "stop_3",
    offset  => "100%",
    'stop-color'  => "rgb(85, 145, 201)",
);

# use explicit element constructor to generate a group element
my $y = $svg->group(
    id => 'group_y',
    style => {
        stroke => 'red',
        fill   => 'black'
    },
);

my $xStart = 10;
my $xScale = 180 / $proteinLen;

my $xPos = $xStart;
my $yPos = 10;
my $yInc = 

$y->rectangle(
    x      => $xPos,
    y      => $yPos,
    width  => 180,
    height => 10,
    rx     => 5.2,
    ry     => 2.4,
    id     => 'protein',
    style => {
        stroke => 'black',
        fill   => 'url(#gradient_1)'
    },

);

my $svgIdx = 1;
foreach my $method ( @methodOrder ) {
  next if ( ! exists $annotData->{$method} );
  my $annots = $annotData->{$method};
  my $color = $colors{$method};
  foreach my $annot ( @{$annots} ) {
    $yPos += 20;
    foreach my $segment ( @{$annot} ) {
      my $segStart = $segment->[0];
      my $segEnd = $segment->[1];
      $xPos = $xStart + ( $segStart * $xScale );
      my $xPos2 = $xStart + ( $segEnd * $xScale );
      my $xLen = $xPos2-$xPos;
      $y->rectangle(
        x      => $xPos,
        y      => $yPos,
        width  => $xLen,
        height => 10,
        #rx     => 5.2,
        #ry     => 2.4,
        id     => "svg_ele_$svgIdx",
      style => {
          stroke => 'grey',
            fill   => $color 
      },
  
    );
    $svgIdx++;
  }
}
}
 
 
# now render the SVG object, implicitly use svg namespace
print $svg->xmlify;

if ( $protName =~ /^(\S+)\#.*$/ ) {
  $protName = $1;
}
open OUT,">$famName-$protName.svg" or die "Could not open $famName-$protName.svg for writing!\n";
print OUT $svg->xmlify;
close OUT;


# or, explicitly use svg namespace and generate a document with its own DTD
#print $svg->xmlify(-namespace=>'svg');
# or, explicitly use svg namespace and generate an inline docunent
#print $svg->xmlify(
#    -namespace => "svg",
#    -pubid => "-//W3C//DTD SVG 1.0//EN",
#    -inline   => 1
#);
}
