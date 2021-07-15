#!/bin/csh

set PROJDIR="paper-data"

# SET Nextflow path
set NEXTFLOW="/usr/local/nextflow/nextflow"
#set NEXTFLOW="/home/rhubley/nextflow-21.04.1/nextflow"

set NEXTFLOW_OPTS=""

# Nextflow run options
set NEXTFLOW_RUN_OPTS=""
#set NEXTFLOW_RUN_OPTS=" --cluster nocona "

set SPECIAL_OPTS=""
#set SPECIAL_OPTS=" -qs 200"

foreach SIM ( DNATransTree-1-Tigger1-R3S-eval DNATransTree-2-Tigger1-R3S-eval )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( DNATransTree-1-Charlie1-R3S-eval DNATransTree-2-Charlie1-R3S-eval )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( LINETree-1-L2-R3S-eval LINETree-2-L2-R3S-eval )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( LINETree-1-CR1-R3S-eval LINETree-2-CR1-R3S-eval )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( DNATransTree-1-Tigger1-R3S-eval  DNATransTree-1-Tigger1-R3S-eval  DNATransTree-1-Tigger1-R3S-eval )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end

foreach SIM ( LINETree-1-L2-R3S-eval  LINETree-1-L2-R3S-eval  LINETree-1-L2-R3S-eval )
  ${NEXTFLOW} ${NEXTFLOW_OPTS} run ${NEXTFLOW_RUN_OPTS} ./generateHTMLViz.nf \
                               --dir ${PROJDIR}/${SIM} \
                               ${SPECIAL_OPTS}
end
