#!/bin/csh

set PROJDIR="paper-data"

rm ${PROJDIR}/cons_stats.txt
rm ${PROJDIR}/sps_stats.txt

foreach SIM ( DNATransTree-1-Tigger1-R3S \
              DNATransTree-2-Tigger1-R3S \
              DNATransTree-2-Charlie1-R3S \
              DNATransTree-1-Charlie1-R3S \
              DNATransTree-1-Tigger1-R3S-gput100-mfl2  DNATransTree-1-Tigger1-R3S-gput1500-mfl2  DNATransTree-1-Tigger1-R3S-gput3000-mfl2 \
              LINETree-1-L2-R3S-gput100-mfl2  LINETree-1-L2-R3S-gput1500-mfl2  LINETree-1-L2-R3S-gput3000-mfl2 \
              LINETree-1-L2-R3S \
              LINETree-2-L2-R3S \
              LINETree-1-CR1-R3S \
              LINETree-2-CR1-R3S )
  #./generateTablesAndGraphs.pl ${PROJDIR}/${SIM}-eval
  ./util/CONS-Significance.py ${PROJDIR}/${SIM}-eval/replicates.csv >> ${PROJDIR}/cons_stats.txt
  ./util/SPS-Significance.py ${PROJDIR}/${SIM}-eval/replicates.csv >> ${PROJDIR}/sps_stats.txt
end



