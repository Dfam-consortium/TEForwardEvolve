#!/bin/csh

set PROJDIR="paper-data"

set NEXTFLOW="/usr/local/nextflow/nextflow"
set NEXTFLOW_RUN_OPTS=""
set NEXTFLOW_OPTS=""

rm -f ${PROJDIR}/naturalSeqsRun.log
#
set EVAL_OPTS="--refiner --muscle --mafft --clustalo --dialign --kalign --fsa --probcons --tcoffee"

foreach SIM ( Arthur2 Tigger10 Zaphod Zaphod2 )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./evalNaturalSeqs.nf \
                               ${EVAL_OPTS} --curated_protein ${PROJDIR}/${SIM}/${SIM}.aa --db_proteins ${PROJDIR}/${SIM}/related_proteins.aa \
                               --familyDir ${PROJDIR}/${SIM} \
                               --outputDir ${PROJDIR}/${SIM}-eval
  cat .nextflow.log >> ${PROJDIR}/naturalSeqsRun.log
end
