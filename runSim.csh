#!/bin/csh

set PROJDIR="paper-data"

###
### Divergence Simulations
###   - Equally distributed fragmentation sizes ranging from 100bp to full length
###   - No minimum full length sequences
###   - Varying "generations per unit time" from 100 to 6000
###
# DNATree-1, Tigger1/Charlie1, Theoretical trevolver matrix
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --transitionFactor 2 --outputDir ${PROJDIR}/DNATransTree-1-Tigger1-K80
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Charlie1.fa --transitionFactor 2 --outputDir ${PROJDIR}/DNATransTree-1-Charlie1-K80
# Tree-2, Tigger1/Charlie1, Theoretical trevolver matrix
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-2.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --transitionFactor 2 --outputDir ${PROJDIR}/DNATransTree-2-Tigger1-K80
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-2.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Charlie1.fa --transitionFactor 2 --outputDir ${PROJDIR}/DNATransTree-2-Charlie1-K80

# DNATree-1, Tigger1/Charlie1, R3S matrix
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --outputDir ${PROJDIR}/DNATransTree-1-Tigger1-R3S
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Charlie1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --outputDir ${PROJDIR}/DNATransTree-1-Charlie1-R3S
# Tree-2, Tigger1/Charlie1, R3S matrix
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-2.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --outputDir ${PROJDIR}/DNATransTree-2-Tigger1-R3S
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-2.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Charlie1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --outputDir ${PROJDIR}/DNATransTree-2-Charlie1-R3S

# LINETree-1, L2/CR1, Theoretical trevolver matrix
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --transitionFactor 2 --outputDir ${PROJDIR}/LINETree-1-L2-K80
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/CR1_Mam.fa --transitionFactor 2 --outputDir ${PROJDIR}/LINETree-1-CR1-K80
# LINETree-2, L2/CR1, Theoretical trevolver matrix
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-2.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --transitionFactor 2 --outputDir ${PROJDIR}/LINETree-2-L2-K80
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-2.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/CR1_Mam.fa --transitionFactor 2 --outputDir ${PROJDIR}/LINETree-2-CR1-K80

# LINETree-1, L2/CR1, R3S matrix
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --outputDir ${PROJDIR}/LINETree-1-L2-R3S
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/CR1_Mam.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --outputDir ${PROJDIR}/LINETree-1-CR1-R3S
# LINETree-2, L2/CR1, R3S matrix
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-2.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --outputDir ${PROJDIR}/LINETree-2-L2-R3S
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-2.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/CR1_Mam.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --outputDir ${PROJDIR}/LINETree-2-CR1-R3S

###
### Fragmentation Simulations
###   - Fragment size determined by a lognormal distribution with min size 50bp and up to full length.
###   - No minimum full length sequences, or minimum 2 full length sequences
###   - Varying fragment size lognormal distribution mean from 75bp to 1200 with a fixed standard deviation of 300
###

# First, no minimum full length sequences
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 100 --outputDir ${PROJDIR}/DNATransTree-1-Tigger1-R3S-gput100

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 1500 --outputDir ${PROJDIR}/DNATransTree-1-Tigger1-R3S-gput1500

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 3000 --outputDir ${PROJDIR}/DNATransTree-1-Tigger1-R3S-gput3000

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 100 --outputDir ${PROJDIR}/LINETree-1-L2-R3S-gput100

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 1500 --outputDir ${PROJDIR}/LINETree-1-L2-R3S-gput1500

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 3000 --outputDir ${PROJDIR}/LINETree-1-L2-R3S-gput3000

# Now, minimum full length of 2 sequences
/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --minFullLen 2 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 100 --outputDir ${PROJDIR}/DNATransTree-1-Tigger1-R3S-gput100-mfl2

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --minFullLen 2 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 1500 --outputDir ${PROJDIR}/DNATransTree-1-Tigger1-R3S-gput1500-mfl2

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --minFullLen 2 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/DNATransTree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/Tigger1.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 3000 --outputDir ${PROJDIR}/DNATransTree-1-Tigger1-R3S-gput3000-mfl2

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --minFullLen 2 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 100 --outputDir ${PROJDIR}/LINETree-1-L2-R3S-gput100-mfl2

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --minFullLen 2 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 1500 --outputDir ${PROJDIR}/LINETree-1-L2-R3S-gput1500-mfl2

/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./genSimulations.nf --minFragLen 50 --minFullLen 2 --tree /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/trees/LINETree-1.nw --seed /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/seeds/L2.fa --matrix /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/matrices/r3s-t1.1.matrix --varyFragSize --gput 3000 --outputDir ${PROJDIR}/LINETree-1-L2-R3S-gput3000-mfl2

