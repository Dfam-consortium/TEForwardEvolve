#!/usr/bin/env nextflow
/*
vim: syntax=groovy

evalMultipleAlignment.nf : Evaluate multiple alignment programs on sequence sets

 Parameters:

     --outputDir <dir>    : Directory to store the results
     --benchmarkDir <dir> : Directory where the benchmark sequences are found
     --seed <file>        : File containing the seed for the simulation
     --muscle             : Run Muscle [optional]
     --refiner            : Run Refiner [optional]
     --clustalw2          : Run ClustalW2 [optional]
     --mafft              : Run Mafft [optiona]
     --cluster            : Either "local", "quanah", "hrothgar" or "griz" 
                            default="local"

 Example:

    nextflow run /blah/evalMultipleAlignment.nf \
                    --benchmarkDir LINETree-1 \
                    --seed L2.fa
                    --muscle \
                    --refiner \
                    --cluster griz

Robert Hubley, 10/2020
*/


// Defaults
params.cluster = "local"
params.muscle = false
params.refiner = false
params.clustalw2 = false
params.mafft = false
params.benchmarkDir = "${workflow.projectDir}/sample"
params.outputDir = "results"
params.seed = "${workflow.projectDir}/sample/L2.fa"

runMuscle = params.muscle
runRefiner = params.refiner
runMafft = params.mafft
runClustalW2 = params.clustalw2
runFastSP = false
runQScore = true
outputDir = params.outputDir

// Default software dependencies ( see localizations in cluster sections )
qscoreDir = "/home/rhubley/projects/DNAMultipleAlignment/qscore"
phrapDir = "/usr/local/phrap"
hmmerDir = "/usr/local/hmmer/bin"
mafftDir = "/usr/local/mafft/bin"
clustalW2Dir = "/usr/local/bin"
dartDir = "/home/rhubley/projects/DNAMultipleAlignment/dart/bin"
muscleDir = "/usr/local/bin"
repeatmodelerDir = "/home/rhubley/projects/RepeatModeler"
fastSPDir = "/home/rhubley/projects/DNAMultipleAlignment/FastSP"


redoFiles = Channel.fromPath( params.outputDir + "/*/gput*-train-*.cons.fa" )


//
// Setup executor for different environments, particularly
// well-known-environments.
//
if ( params.cluster == "local" ) {
  thisExecutor = "local"
  thisQueue = ""
  thisOptions = ""
  thisScratch = false
}else if ( params.cluster == "quanah" || params.cluster == "hrothgar" ){
  proc = 12
  thisExecutor = "sge"
  thisQueue = "omni"
  thisOptions = "-pe sm ${proc} -P quanah -S /bin/bash"
  ucscToolsDir="/home/rhubley/ucscTools"
  repeatMaskerDir="/home/rhubley/RepeatMasker"
  thisScratch = false
}else if ( params.cluster == "griz" ) {
  proc = 12
  thisExecutor = "slurm"
  //thisQueue = "wheeler_lab_small_cpu"
  thisQueue = "wheeler_lab_large_cpu"
  thisOptions = "--tasks=1 --cpus-per-task=${proc}"
  ucscToolsDir="/home/rh105648e/ucscTools"
  repeatMaskerDir="/home/rh105648e/RepeatMasker-4.1.1"
  thisScratch = "/state/partition1"
}

log.info "evalMultipleAlign.nf : Multiple Alignment Evaluation ver 0.1"
log.info "===================================================================="
log.info "working directory   : " + workflow.workDir
log.info "RepeatModeler DIR   : " + repeatmodelerDir
log.info "Benchmark DIR       : " + params.benchmarkDir
log.info "Seed File           : " + params.seed
log.info "Muscle DIR          : " + muscleDir
log.info "Mafft DIR           : " + mafftDir
log.info "ClustalW2 DIR       : " + clustalW2Dir
log.info "Cluster             : " + params.cluster
log.info "Output DIR          : " + outputDir
log.info "Queue/Partititon    : " + thisQueue



process doitprocess {
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  file(consFile) from redoFiles

  output:
  file "*nw_score"
  

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (consFile.name =~ /(gput\d+-train)/)[0][0]
  gputPrefix = (consFile.name =~ /(gput\d+)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = consFile.toRealPath().getName(consFile.toRealPath().getNameCount() - 2)
  """
  ${workflow.projectDir}/compareToTestNeedle.pl ${consFile} /home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/${params.benchmarkDir}/${repDir}/${gputPrefix}-test-seqs.fa > ${consFile.baseName}.nw_score
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


