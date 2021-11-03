
# TEForwardEvolve

## Transposable Element Sequence Simulation

This project contains perl/python/Nextflow code for a foward evolution DNA
sequence simulator.  This was used to develop a Transposable Element benchmark 
for evaluating multiple sequence alignment tools.  Please see the preprint
paper at:

https://www.biorxiv.org/content/10.1101/2021.08.17.456740v1


### Dependencies

* Nextflow 21.04.1 or higher
* Python 3
  * SciPy module
  * Numpy module
  * Matplotlib module
* Perl 5.16 or higher
  * SVG module
     
### Examples

Simple parameters:
```
TEForwardEvolve.pl -matrix matrices/r3s-t1.1.matrix  \
                   -tree trees/DNATransTree-1.nw \
                   -seed seeds/Tigger1.fa \
                   -generations_per_unit_time 1750 \
                   -verbosity 2
```

With Fragmentation Parameters:
```
TEForwardEvolve.pl -matrix matrices/r3s-t1.1.matrix  \
                   -tree trees/DNATransTree-1.nw \
                   -seed seeds/Tigger1.fa \
                   -fragment_size_mean 75 \ 
                   -fragment_size_stdev 300 \
                   -min_full_len 2 \
                   -min_frag_len 50 \
                   -generations_per_unit_time 1750 \
                   -verbosity 2
```

### Files/Directories

```
seeds/                A folder containing sequences used to seed simulations
  L2.fa                   LINE L2 family
  CR1_Mam.fa              LINE CR1 family
  Tigger1.fa              DNA Transposon Tigger1 family
  Charlie1.fa             DNA Transposon Charlie1 family
trees/                Simulated trees generated by util/genRandomTETrees.pl
  DNATransTree-1.nw       DNA Transposon tree simulation #1 used in the main text
  DNATransTree-2.nw       DNA Transposon tree simulation #2 used in the supplemental
  LINETree-1.nw           LINE tree simulation #1 used in the main text
  LINETree-2.nw           LINE tree simulation #2 used in the supplelmental
matrices/             Substitution matrices used in sequence simulations
  r3s-t1.1.mod            Instantaneous rate matrix 64x64 for mammals supplied by Adam Siepel
  r3s-t1.1.matrix         Instantaneous rate matrix 64x4 with fixed flanking bases derived
                            from r3s-1.1.mod

TEForwardEvolve.pl    The sequence evolution main program

runSim.csh            Script to generate the collection of simulated sequences and reference
                         MSAs.
runEval.csh           Script to run all MSA tools on the simulated sequences and run the
                         various evaluation programs.
runStat.csh           Script to collect evaluation results into graphs, tables, and stats.
```

The remaining files are called by the above scripts.

