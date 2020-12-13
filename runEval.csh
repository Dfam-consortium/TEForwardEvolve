#!/bin/csh

set PROJDIR="paper-data"

# Important to keep the work directories on a local disk.
# On Dfam-Dev:
set NEXTFLOW_OPTS="-log /local/MSA/my.log"
set NEXTFLOW_RUN_OPTS="-w /local/MSA/work"

#
set EVAL_OPTS="--refiner --muscle --mafft --clustalw2 --dialign --kalign"


#foreach SIM ( DNATransTree-1-Tigger1-K80  DNATransTree-1-Tigger1-R3S \
#              DNATransTree-2-Tigger1-K80  DNATransTree-2-Tigger1-R3S )
#  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
#                               ${EVAL_OPTS} --seed seeds/Tigger1.fa \
#                               --benchmarkDir ${PROJDIR}/${SIM} \
#                               --outputDir ${PROJDIR}/${SIM}-eval
#end
#
#foreach SIM ( DNATransTree-2-Charlie1-K80 DNATransTree-2-Charlie1-R3S \
#              DNATransTree-1-Charlie1-K80 DNATransTree-1-Charlie1-R3S )
#  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
#                               ${EVAL_OPTS} --seed seeds/Charlie1.fa \
#                               --benchmarkDir ${PROJDIR}/${SIM} \
#                               --outputDir ${PROJDIR}/${SIM}-eval
#end

foreach SIM ( LINETree-1-L2-K80  LINETree-1-L2-R3S \
              LINETree-2-L2-K80  LINETree-2-L2-R3S )
  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
                               ${EVAL_OPTS} --seed seeds/L2.fa \
                               --benchmarkDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval
end

foreach SIM ( LINETree-1-CR1-K80 LINETree-1-CR1-R3S \
              LINETree-2-CR1-K80 LINETree-2-CR1-R3S )
  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
                               ${EVAL_OPTS} --seed seeds/CR1_Mam.fa \
                               --benchmarkDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval

end

