#!/usr/bin/env nextflow
/*
vim: syntax=groovy

genSimulations.nf : Run TEFowardEvolve with a single tree in replicate

 Parameters:

     --outputDir <dir>      : Directory to store results in
     --tree <tree_file>     : Newick file for the simulation
     --seed <seed_file>     : Seed sequence FASTA file
     --matrix <matrix_file> : Matrix
     --minFragLen <#>       : Minimum fragment length
     --minFullLen <#>       : Minimu number of full length seqs
     --transitionFactor <#> : transition factor
     --cpgFactor <#>        : CpG factor
     --varyFragSize         : Vary fragment size rather than gput scales [ default: off ]
     --gput <#>             : If varying fragment sizes this sets the fixed gput [default: 100]
     --replicates <#>       : Number of replicates to generate [default 10]
     --cluster <name>       : Either "local", "quanah", "hrothgar" or "griz"

 Example:

    nextflow run /blah/genSimulations.nf \
                    --outputDir LINETree-1 \
                    --cluster griz

Robert Hubley, 10/2020
*/


// Defaults
params.cluster = "local"
params.seed = "${workflow.projectDir}/seeds/L2.fa"
params.tree = "${workflow.projectDir}/trees/LINETree.nw"
params.minFullLen = "undefined"
params.matrix = "undefined"
params.minFragLen = "undefined"
params.transitionFactor = "undefined"
params.cpgFactor = "undefined"
params.replicates = 10
params.gput = 100
params.varyFragSize = false
params.outputDir = "${workflow.projectDir}/ExampleSim-LINETree-L2"

treeValue = Channel.value(params.tree)
seedValue = Channel.value(params.seed)

matrixParam = ""
if ( params.matrix != "undefined" ) {
  matrixParam = "-matrix ${params.matrix}"
}else{
  if ( params.transitionFactor != "undefined" ) {
    matrixParam = "-transition_factor ${params.transitionFactor}"
  }
  if ( params.cpgFactor != "undefined" ) {
    matrixParam = matrixParam + " -cpg_factor ${params.cpgFactor}"
  }
}

otherParam = ""
if ( params.minFullLen != "undefined" ) {
  otherParam += " -min_full_len ${params.minFullLen}"
}
if ( params.minFragLen != "undefined" ) {
  otherParam += " -min_frag_len ${params.minFragLen}"
}

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
  thisScratch = false
}else if ( params.cluster == "griz" ) {
  proc = 12
  thisExecutor = "slurm"
  //thisQueue = "wheeler_lab_small_cpu"
  thisQueue = "wheeler_lab_large_cpu"
  thisOptions = "--tasks=1 --cpus-per-task=${proc}"
  thisScratch = "/state/partition1"
}

log.info "genSimulations.nf : Generate replicates for a single simulation tree ver 0.1"
log.info "============================================================================="
log.info "Working Directory   : " + workflow.workDir
log.info "Newick Tree         : " + params.tree
log.info "Seed File           : " + params.seed
if ( matrixParam == "" ) {
  log.info "Matrix Param        : Default TRevolver Matrix"
}else {
  log.info "Matrix Param        : " + matrixParam
}
if ( otherParam != "" ) {
  log.info "Other Param         : " + otherParam
}
log.info "Replicates          : " + params.replicates
if ( params.varyFragSize ) {
log.info "Parameter           : Fragmentation Size Distribution [gput=" + params.gput + "]"
}else {
log.info "Parameter           : Generations Per Unit Time"
}
log.info "Output Directory    : " + params.outputDir
log.info "Cluster             : " + params.cluster
log.info "Queue/Partititon    : " + thisQueue

outputDir = params.outputDir

setChan = Channel.from( "test", "train" )
replicateChan = Channel.of(1..params.replicates)

scaleSimulationChan = Channel.from()
fragSimulationChan = Channel.from()
if ( params.varyFragSize ) {
  fragChan = Channel.from([75, 300], [100, 300], [150,300], [200,300], [250,300], [300, 300], [400, 300], [500,300], [600,300], [800,300], [1000,300], [1200, 300])
  fragSimulationChan = setChan.combine(fragChan.combine(replicateChan))
}else {
  scaleChan = Channel.from(100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 4000, 5000, 6000) 
  scaleSimulationChan = setChan.combine(scaleChan.combine(replicateChan))
}

//setChan = Channel.from( "test", "train" )
//scaleChan = Channel.from(100, 250)
//replicateChan = Channel.of(1..2)




process runTEForwardEvolveScales {
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "rep-${replicate}/$filename" }

  input:
  path seedFile from seedValue
  path treeFile from treeValue
  set set_val, scale, replicate from scaleSimulationChan

  output:
  file "*.log"
  file "*refmsa.fa"
  file "*seqs.fa"

  script:
  //log.info "Running: " + seedFile.name + " set= " + set_val + " scale=" + scale + " replicate=" + replicate
  """
  ${workflow.projectDir}/TEForwardEvolve.pl ${matrixParam} ${otherParam} -tree ${treeFile} -seed ${seedFile} -generations_per_unit_time ${scale} -verbosity 2 >& out.log
  mv out.log gput${scale}-${set_val}.log
  mv output-msa.fa gput${scale}-${set_val}-refmsa.fa
  mv output-seqs.fa gput${scale}-${set_val}-seqs.fa
  """
}

process runTEForwardEvolveFrags {
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "rep-${replicate}/$filename" }

  input:
  path seedFile from seedValue
  path treeFile from treeValue
  set set_val, fmean, fstdev, replicate from fragSimulationChan

  output:
  file "*.log"
  file "*refmsa.fa"
  file "*seqs.fa"

  script:
  //log.info "Running: " + seedFile.name + " set= " + set_val + " fmean=" + fmean + " fstdev=" + fstdev + " replicate=" + replicate
  """
  ${workflow.projectDir}/TEForwardEvolve.pl ${matrixParam} ${otherParam} -fragment_size_mean ${fmean} -fragment_size_stdev ${fstdev} -tree ${treeFile} -seed ${seedFile} -generations_per_unit_time ${params.gput} -verbosity 2 >& out.log
  mv out.log frag${fmean}-${fstdev}-${set_val}.log
  mv output-msa.fa frag${fmean}-${fstdev}-${set_val}-refmsa.fa
  mv output-seqs.fa frag${fmean}-${fstdev}-${set_val}-seqs.fa
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

