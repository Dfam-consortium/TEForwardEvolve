#!/usr/bin/env nextflow
/*
vim: syntax=groovy

generateHTMLViz : Generate visualization for MSAs

 Parameters:

     --dir <dir>    : Directory to find input and store output
     --cluster      : Either "local", "quanah", "hrothgar" or "griz" 
                        default="local"

 Example:

    nextflow run /blah/generateHTMLViz.nf \
                    --dir paper-data/DNATransTree-1-Tigger1-R3S-eval
                    --cluster quanah

Robert Hubley, 7/2021
*/


// Defaults
params.cluster = "local"
params.dir = "${workflow.projectDir}/sample"
dataDir = params.dir

// Viz software through RepeatModeler package [ https://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.2a.tar.gz ]
repeatmodelerDir = "/home/rhubley/RepeatModeler-2.0.2a"
//repeatmodelerDir = "/usr/local/RepeatModeler-2.0.2a"

//
// Setup executor for different environments, particularly
// well-known-environments.
//
if ( params.cluster == "local" ) {
  thisExecutor = "local"
  thisQueue = ""
  thisOptions = ""
  thisScratch = false
}else if ( params.cluster == "nocona" ){
  // nocona:
  //   240 nodes, 38720 cores
  //   128 cores/node, 512 GB/node
  //   4 GB/core
  thisExecutor = "slurm"
  thisQueue = "nocona"
  thisOptions = "--ntasks=1 --nodes=1"
  thisScratch = false
}else if ( params.cluster == "griz" ) {
  proc = 12
  thisExecutor = "slurm"
  thisQueue = "wheeler_lab_large_cpu"
  thisOptions = "--tasks=1 --cpus-per-task=${proc}"
  thisScratch = "/state/partition1"
}

log.info "generateHTMLViz.nf : MSA Visualization ver 0.1"
log.info "==============================================="
log.info "working directory   : " + workflow.workDir
log.info "RepeatModeler DIR   : " + repeatmodelerDir
log.info "Data DIR            : " + dataDir
log.info "Cluster             : " + params.cluster
log.info "Queue/Partititon    : " + thisQueue

//workChan = Channel.fromPath("paper-data/test-eval/rep-*/*-train-{muscle,mafft,refiner,clustalo,dialign,kalign,fsa,refmsa}.fa").view()
workChan = Channel.fromPath("${dataDir}/rep-*/*-train-{muscle,mafft,refiner,clustalo,dialign,kalign,fsa,refmsa}.fa")
       .map { it -> tokens = (it =~ /(rep-\d+)\/((gput|frag)\d+)/)[0]; [tokens[1], tokens[2], it ] }
       .groupTuple(by:[0,1], sort: true)


process visualizeMSA{
  cpus 1
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch
  publishDir "${dataDir}/html", mode: 'copy'

  input:
  set rep, param, file(MSAFiles) from workChan
  
  output: 
  file "rep*.html"

  script:
  """
  ${repeatmodelerDir}/util/viewMultipleMSA.pl ${MSAFiles}
  mv MultMSA.html rep-${rep}-${param}.html
  ${repeatmodelerDir}/util/viewMultipleMSA.pl -fullmsa ${MSAFiles}
  mv MultMSA.html rep-${rep}-${param}-fullmsa.html
  """
}

workflow.onComplete {
            log.info """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
}
   

