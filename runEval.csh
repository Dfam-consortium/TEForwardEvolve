#!/bin/csh

set PROJDIR="paper-data"

# Important to keep the work directories on a local disk.
# On Dfam-Dev:
set NEXTFLOW_OPTS="-log /local/MSA/my.log"
set NEXTFLOW_RUN_OPTS="-w /local/MSA/work"

#
set EVAL_OPTS="--refiner --muscle --mafft --clustalo --dialign --kalign --fsa"

foreach SIM ( DNATransTree-1-Tigger1-R3S DNATransTree-2-Tigger1-R3S )
  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
                               ${EVAL_OPTS} --seed seeds/Tigger1.fa \
                               --benchmarkDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval
end

foreach SIM ( DNATransTree-2-Charlie1-R3S DNATransTree-1-Charlie1-R3S )
  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
                               ${EVAL_OPTS} --seed seeds/Charlie1.fa \
                               --benchmarkDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval
end

foreach SIM ( LINETree-1-L2-R3S LINETree-2-L2-R3S )
  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
                               ${EVAL_OPTS} --seed seeds/L2.fa \
                               --benchmarkDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval
end

foreach SIM ( LINETree-1-CR1-R3S LINETree-2-CR1-R3S )
  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
                               ${EVAL_OPTS} --seed seeds/CR1_Mam.fa \
                               --benchmarkDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval
end

foreach SIM ( DNATransTree-1-Tigger1-R3S-gput100-mfl2  DNATransTree-1-Tigger1-R3S-gput1500-mfl2  DNATransTree-1-Tigger1-R3S-gput3000-mfl2 )
  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
                               ${EVAL_OPTS} --seed seeds/Tigger1.fa \
                               --benchmarkDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval
end

foreach SIM ( LINETree-1-L2-R3S-gput100-mfl2  LINETree-1-L2-R3S-gput1500-mfl2  LINETree-1-L2-R3S-gput3000-mfl2 )
  /usr/local/nextflow/nextflow ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalMultipleAlign.nf \
                               ${EVAL_OPTS} --seed seeds/L2.fa \
                               --benchmarkDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval
end

