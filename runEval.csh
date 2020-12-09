#!/bin/csh

#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Tigger1.fa \
#  --benchmarkDir run2/DNATransTree-1-Tigger1-K80 \
#  --outputDir run2/DNATransTree-1-Tigger1-K80-eval
#
#/usr/local/nextflow/nextflow run ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Charlie1.fa \
#  --benchmarkDir run2/DNATransTree-1-Charlie1-K80 \
#  --outputDir run2/DNATransTree-1-Charlie1-K80-eval
#
#/usr/local/nextflow/nextflow run ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Tigger1.fa \
#  --benchmarkDir run2/DNATransTree-2-Tigger1-K80 \
#  --outputDir run2/DNATransTree-2-Tigger1-K80-eval
#
#/usr/local/nextflow/nextflow run ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Charlie1.fa \
#  --benchmarkDir run2/DNATransTree-2-Charlie1-K80 \
#  --outputDir run2/DNATransTree-2-Charlie1-K80-eval
#
#/usr/local/nextflow/nextflow run ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Tigger1.fa \
#  --benchmarkDir run2/DNATransTree-1-Tigger1-R3S \
#  --outputDir run2/DNATransTree-1-Tigger1-R3S-eval

#/usr/local/nextflow/nextflow run ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Tigger1.fa \
#  --benchmarkDir run2/DNATransTree-1-Tigger1-R3S \
#  --outputDir run2/DNATransTree-1-Tigger1-R3S-evalblk

/usr/local/nextflow/nextflow run ./evalMultipleAlign.nf \
  --dialign \
  --kalign \
  --seed seeds/Tigger1.fa \
  --benchmarkDir run2/DNATransTree-1-Tigger1-R3S \
  --outputDir run2/DNATransTree-1-Tigger1-R3S-added

#
#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Tigger1.fa \
#  --benchmarkDir run2/DNATransTree-2-Tigger1-R3S \
#  --outputDir run2/DNATransTree-2-Tigger1-R3S-eval
#
#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Charlie1.fa \
#  --benchmarkDir run2/DNATransTree-1-Charlie1-R3S \
#  --outputDir run2/DNATransTree-1-Charlie1-R3S-eval
#
#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/Charlie1.fa \
#  --benchmarkDir run2/DNATransTree-2-Charlie1-R3S \
#  --outputDir run2/DNATransTree-2-Charlie1-R3S-eval
#
#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/L2.fa \
#  --benchmarkDir run2/LINETree-1-L2-K80 \
#  --outputDir run2/LINETree-1-L2-K80-eval
#
#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/L2.fa \
#  --benchmarkDir run2/LINETree-1-L2-R3S \
#  --outputDir run2/LINETree-1-L2-R3S-eval
#
#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/L2.fa \
#  --benchmarkDir run2/LINETree-2-L2-K80 \
#  --outputDir run2/LINETree-2-L2-K80-eval
#
#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/L2.fa \
#  --benchmarkDir run2/LINETree-2-L2-R3S \
#  --outputDir run2/LINETree-2-L2-R3S-eval

#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/CR1.fa \
#  --benchmarkDir run2/LINETree-1-CR1-K80 \
#  --outputDir run2/LINETree-1-CR1-K80-eval

#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/CR1.fa \
#  --benchmarkDir run2/LINETree-1-CR1-R3S \
#  --outputDir run2/LINETree-1-CR1-R3S-eval

#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/CR1.fa \
#  --benchmarkDir run2/LINETree-2-CR1-K80 \
#  --outputDir run2/LINETree-2-CR1-K80-eval

#/usr/local/nextflow/nextflow -log /local/MSA/my.log run -w /local/MSA/work ./evalMultipleAlign.nf \
#  --refiner \
#  --muscle \
#  --mafft \
#  --clustalw2 \
#  --seed seeds/CR1.fa \
#  --benchmarkDir run2/LINETree-2-CR1-R3S \
#  --outputDir run2/LINETree-2-CR1-R3S-eval

