#!/bin/csh

set PROJDIR="paper-data"

foreach SIM ( DNATransTree-1-Tigger1-K80  DNATransTree-1-Tigger1-R3S \
              DNATransTree-2-Tigger1-K80  DNATransTree-2-Tigger1-R3S \
              DNATransTree-2-Charlie1-K80 DNATransTree-2-Charlie1-R3S \
              DNATransTree-1-Charlie1-K80 DNATransTree-1-Charlie1-R3S \
              DNATransTree-1-Tigger1-R3S-gput100  DNATransTree-1-Tigger1-R3S-gput1500  DNATransTree-1-Tigger1-R3S-gput3000 \
              DNATransTree-1-Tigger1-R3S-gput100-mfl2  DNATransTree-1-Tigger1-R3S-gput1500-mfl2  DNATransTree-1-Tigger1-R3S-gput3000-mfl2 \
              LINETree-1-L2-R3S-gput100  LINETree-1-L2-R3S-gput1500  LINETree-1-L2-R3S-gput3000 \
              LINETree-1-L2-R3S-gput100-mfl2  LINETree-1-L2-R3S-gput1500-mfl2  LINETree-1-L2-R3S-gput3000-mfl2 \
              LINETree-1-L2-K80  LINETree-1-L2-R3S \
              LINETree-2-L2-K80  LINETree-2-L2-R3S \
              LINETree-1-CR1-K80 LINETree-1-CR1-R3S \
              LINETree-2-CR1-K80 LINETree-2-CR1-R3S )
  ./generateTablesAndGraphs.pl ${PROJDIR}/${SIM}-eval
end



