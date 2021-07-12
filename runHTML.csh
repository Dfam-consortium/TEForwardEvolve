#!/bin/csh

set PROJDIR="paper-data"

# SET Nextflow path
#set NEXTFLOW="/usr/local/nextflow/nextflow"
set NEXTFLOW="/home/rhubley/nextflow-21.04.1/nextflow"

set NEXTFLOW_OPTS=""

# Nextflow run options
#set NEXTFLOW_RUN_OPTS=""
set NEXTFLOW_RUN_OPTS=" --cluster nocona "

set SPECIAL_OPTS=" -qs 200"

#
set EVAL_OPTS="--refiner --muscle --mafft --clustalo --dialign --kalign --fsa"

foreach SIM ( DNATransTree-1-Tigger1-R3S DNATransTree-2-Tigger1-R3S )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( DNATransTree-1-Charlie1-R3S DNATransTree-2-Charlie1-R3S )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( LINETree-1-L2-R3S LINETree-2-L2-R3S )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( LINETree-1-CR1-R3S LINETree-2-CR1-R3S )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( DNATransTree-1-Tigger1-R3S-gput100-mfl2  DNATransTree-1-Tigger1-R3S-gput1500-mfl2  DNATransTree-1-Tigger1-R3S-gput3000-mfl2 )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( LINETree-1-L2-R3S-gput100-mfl2  LINETree-1-L2-R3S-gput1500-mfl2  LINETree-1-L2-R3S-gput3000-mfl2 )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end
