#!/bin/csh

set PROJDIR="paper-data"

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