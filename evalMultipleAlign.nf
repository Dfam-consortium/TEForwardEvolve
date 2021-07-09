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
     --dialign            : Run DialignTX [optional]
     --kalign             : Run Kalign [optional]
     --fsa                : Run FSA [optional]
     --clustalw2          : Run ClustalW2 [optional]
     --clustalo           : Run ClustalOmega [optional]
     --opal               : Run Opal [optional]
     --mafft              : Run Mafft [optional]
     --cmConsEval         : Run crossmatch consensus evaluation [optional]
     --nhmmerConsEval     : Run nhmmer consensus evaluation [optional]
     --nhmmerHMMEval      : Run nhmmer HMM evaluation [optional]
     --dartScore          : Run additional DART cmpalign score [optional]
     --cluster            : Either "local", "quanah", "hrothgar" or "griz" 
                            default="local"

 Example:

    nextflow run /blah/evalMultipleAlignment.nf \
                    --benchmarkDir LINETree-1 \
                    --seed L2.fa \
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
params.fsa = false
params.dialign = false
params.kalign = false
params.clustalo = false
params.opal = false
params.cmConsEval = false
params.nhmmerConsEval = false
params.nhmmerHMMEval = false
params.dartScore = false
params.benchmarkDir = "${workflow.projectDir}/sample"
params.outputDir = "results"
params.seed = "${workflow.projectDir}/sample/L2.fa"

runMuscle = params.muscle
runRefiner = params.refiner
runMafft = params.mafft
runClustalW2 = params.clustalw2
runDialign = params.dialign
runKalign = params.kalign
runFSA = params.fsa
runClustalO = params.clustalo
runOpal = params.opal
runFastSP = false
runQScore = true
runCMConsEval = params.cmConsEval
runNhmmerConsEval = params.nhmmerConsEval
runNhmmerHMMEval = params.nhmmerHMMEval
runDartScore = params.dartScore
outputDir = params.outputDir
// This seems to be important for this special variable
// othewise the hash'es for the resume do not match. It's 
// both necessary to have a new variable rather than use
// ${workflow.projectDir} in a process definition *and*
// "def" is necessary.
def projDir = workflow.projectDir

// Default software dependencies ( see localizations in cluster sections )

// SPS calculation through qscore program
qscoreDir = "/home/rhubley/bin"
//qscoreDir = "/home/rhubley/projects/DNAMultipleAlignment/qscore"

// SPS calcuation throuh dart AMA program
dartDir = "/home/rhubley/bin"
//dartDir = "/home/rhubley/projects/DNAMultipleAlignment/dart/bin"

// SPS calcuation through fastSP package ( removed )
//fastSPDir = "/home/rhubley/projects/DNAMultipleAlignment/FastSP"

// Phil Greens crossmatch program - for consensus model eval
phrapDir = "/home/rhubley/phrap-1.090518"
//phrapDir = "/usr/local/phrap"

// Eddy/Wheeler nhmmer program - for HMM model eval
hmmerDir = "/home/rhubley/hmmer-3.3.2/bin"
//hmmerDir = "/usr/local/hmmer/bin"

// 
// Aligners
//
// MAFFT aligner [ https://mafft.cbrc.jp/alignment/software/mafft-7.481-without-extensions-src.tgz ]
//    vi core/Makefile -- set prefix correctly otherwise hardcoded paths don't get set correctly
//    cd core; make; make install
mafftDir = "/home/rhubley/mafft-7.481/bin"
//mafftDir = "/usr/local/mafft/bin"

// DIALIGN aligner [ http://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.tar.gz ]
//    CPPFLAGS=-O3 -funroll-loops  -mfpmath=sse -msse  -mmmx
//    NOTE: This is the directory containing the subdirctory bin/ and conf/
dialignDir = "/home/rhubley/DIALIGN-TX_1.0.2"
//dialignDir = "/usr/local/dialign-tx-1.0.2"

// Kalign Aligner 2.0.4 [ http://msa.sbc.su.se/downloads/kalign/current.tar.gz ]
//    mkdir kalign-2.0.4; cd kalign-2.0.4; tar zxvf ../current.tar.gz
kalignDir = "/home/rhubley/kalign-2.0.4"
//kalignDir = "/usr/local/kalign2"

// --unused--
clustalW2Dir = "/usr/local/bin"

// Clustal Omega 1.2.4 [ http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64 ]
//    mkdir clustal-omega-1.2.4; mv clustalo-1.2.4-Ubuntu-x86_64 clustal-omega-1.2.4/clustalo
//    chmod 755 clustal-omega-1.2.4/clustalo
clustalOmegaDir = "/home/rhubley/clustal-omega-1.2.4"
//clustalOmegaDir = "/u1/local/clustal-omega-1.2.4-binary"

// --unused-- Opal Aligner 2.1.3 [ http://opal.cs.arizona.edu/old_distros/opal_2.1.3.tgz ]
opalDir = "/u1/local/opal_2.1.3"

// FSA Aligner 1.15.9 [ https://sourceforge.net/projects/fsa/files/fsa-1.15.9.tar.gz/download ]
//     Depends on MUMmer [ https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download ]
//       mv download MUMmer-3.23.tar.gz; tar zxvf MUMmer-3.23.tar.gz; cd MUMmer-3.23/; make install
//     mv download fsa-1.15.9.tar.gz; tar zxvf fsa-1.15.9.tar.gz; cd fsa-1.15.9; ./configure --with....; make
//     *binaries are in fsa-1.15.9/src/main*
fsaDir = "/home/rhubley/fsa-1.15.9/bin"
//fsaDir = "/usr/local/fsa/bin"

// Muscle 3.8.31 [ https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz ]
//     tar zxvf muscle3.8.31_i86linux64.tar.gz; mkdir muscle-3.8.31; mv muscle3.8.31_i86linux64 muscle-3.8.31/muscle
//     chmod 755 muscle-3.8.31/muscle
muscleDir = "/home/rhubley/muscle-3.8.31"
//muscleDir = "/usr/local/bin"

// Refiner alignment through RepeatModeler package [ https://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.2a.tar.gz ]
repeatmodelerDir = "/home/rhubley/RepeatModeler-2.0.2a"
//repeatmodelerDir = "/usr/local/RepeatModeler-2.0.2a"

//FOR DEBUGGING...limit the files run
//Channel.fromFilePairs( params.benchmarkDir + "/rep-1/gput100-{train-seqs,train-refmsa,test-seqs}.fa", size: 3, flat:true )
//     .into { benchmarkFilesForComp }

//
// E.g "gput100", "gput100-train-seqs.fa", "gput100-train-refmsa.fa", "gput100-test-seqs.fa"
Channel.fromFilePairs( params.benchmarkDir + "/*/*-{train-seqs,train-refmsa,test-seqs}.fa", size: 3, flat:true )
       .into { benchmarkFilesForComp }

// Constant seedFile for use in all processes
seedFile = file(params.seed)

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

log.info "evalMultipleAlign.nf : Multiple Alignment Evaluation ver 0.3"
log.info "============================================================"
log.info "working directory   : " + workflow.workDir
log.info "RepeatModeler DIR   : " + repeatmodelerDir
log.info "Benchmark DIR       : " + params.benchmarkDir
log.info "Seed File           : " + params.seed
if ( runMuscle ) {
  log.info "Muscle DIR          : " + muscleDir
}
if ( runMafft ) {
  log.info "Mafft DIR           : " + mafftDir
}
if ( runClustalW2 ) {
  log.info "ClustalW2 DIR       : " + clustalW2Dir
}
if ( runRefiner ) {
  log.info "Refiner DIR         : " + repeatmodelerDir
}
if ( runDialign ) {
  log.info "Dialign DIR         : " + dialignDir
}
if ( runKalign ) {
  log.info "Kalign DIR          : " + kalignDir
}
if ( runFSA ) {
  log.info "FSA DIR             : " + fsaDir
}
if ( runClustalO ) { 
  log.info "ClustalOmega DIR    : " + clustalOmegaDir
}
if ( runOpal ) {
  log.info "Opal DIR            : " + opalDir
}
log.info "Cluster             : " + params.cluster
log.info "Output DIR          : " + outputDir
log.info "Queue/Partititon    : " + thisQueue


process processRefMSA {
  cpus 2
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  file seedFile
  set gput, file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForComp

  output:
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForClustalW2 
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForMafft
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForMuscle
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForRefiner
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForKalign
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForDialign
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForFSA
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForOpal
  tuple file("*cons.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into benchmarkFilesForClustalO
  file "*cons.fa"
  file "*cons.vs_self"
  file "*cons.vs_seed"
  file "*avgKDiv"
  
  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  ###########
  # Generate consensus from reference MSA to have an practical maximum acheivable given this simulation.
  #   E.g. Linup -consensus -name gput100-train-refmsa gput100-train-refmsa.fa > gput100-train-refmsa.cons.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${referenceMSAFile} > ${referenceMSAFile.baseName}.cons.fa
  ${repeatmodelerDir}/util/Linup ${referenceMSAFile} | grep "Avg Kimura Div" > ${referenceMSAFile.baseName}.avgKDiv
  # Compare consensus to input SEED sequence and evaluate
  #   E.g. compareConsensiNeedle.pl gput100-train-refmsa.cons.fa L2.fa > gput100-train-refmsa.cons.vs_seed
  ${projDir}/util/compareConsensiNeedle.pl ${referenceMSAFile.baseName}.cons.fa ${seedFile} > ${referenceMSAFile.baseName}.cons.vs_seed
  ${projDir}/util/compareConsensiNeedle.pl ${referenceMSAFile.baseName}.cons.fa ${referenceMSAFile.baseName}.cons.fa > ${referenceMSAFile.baseName}.cons.vs_self
  ###########
  """
}


process runFSA {
  cpus 2
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForFSA

  when:
  runFSA

  output:
  tuple file("*fsa.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into fsaToAMAChan
  tuple file("*fsa.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into fsaToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into fsaToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into fsaToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) into fsaToNhmmerHMMChan mode flatten
  file "*-fsa.fa"
  file "*-fsa.hmm"
  file "*-fsa.cons.fa"
  file "*-fsa.trimmed-cons.fa"
  file "*-fsa.trimmed.hmm"
  file "*-fsa.cons.vs_refmsacons"
  file "*-fsa.trimmed-cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run FSA
  ${fsaDir}/fsa ${referenceSeqFile} > ${simPrefix}-fsa.fa

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-fsa.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-fsa.fa > ${simPrefix}-fsa.cons.fa
  ${projDir}/util/trimUnalignedEdges.pl ${simPrefix}-fsa.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-fsa.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-fsa.cons.fa L2.fa > gput100-train-fsa.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-fsa.trimmed-cons.fa L2.fa > gput100-train-fsa.trimmed-cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-fsa.cons.fa ${referenceMSACons} > ${simPrefix}-fsa.cons.vs_refmsacons 
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-fsa.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-fsa.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-fsa.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-fsa.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-fsa.hmm normal.stk > hmmbuild.log
  """
}

//  KAlign: https://msa.sbc.su.se/cgi-bin/msa.cgi
//    parameters suggested by website for DNA
process runKalign {
  cpus 2  
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForKalign

  when:
  runKalign

  output:
  tuple file("*kalign.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into kalignToAMAChan
  tuple file("*kalign.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into kalignToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into kalignToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into kalignToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) into kalignToNhmmerHMMChan mode flatten
  file "*-kalign.fa"
  file "*-kalign.hmm"
  file "*-kalign.cons.fa"
  file "*-kalign.trimmed-cons.fa"
  file "*-kalign.trimmed.hmm"
  file "*-kalign.cons.vs_refmsacons"
  file "*-kalign.trimmed-cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run Kalign
  ${kalignDir}/kalign -gpo 80 -gpe 3 -tgpe 3 -bonus 0 ${referenceSeqFile} > ${simPrefix}-kalign.fa

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-kalign.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-kalign.fa > ${simPrefix}-kalign.cons.fa
  ${projDir}/util/trimUnalignedEdges.pl ${simPrefix}-kalign.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-kalign.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-kalign.cons.fa L2.fa > gput100-train-kalign.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-kalign.trimmed-cons.fa L2.fa > gput100-train-kalign.trimmed-cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-kalign.cons.fa ${referenceMSACons} > ${simPrefix}-kalign.cons.vs_refmsacons 
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-kalign.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-kalign.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-kalign.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-kalign.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-kalign.hmm normal.stk > hmmbuild.log
  """
}

//  DIALIGN: http://dialign.gobics.de/
// Switched to dialigntx because it does appear to perform better
process runDialign {
  cpus 2
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForDialign

  when:
  runDialign

  output:
  tuple file("*dialign.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into dialignToAMAChan
  tuple file("*dialign.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into dialignToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into dialignToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into dialignToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) into dialignToNhmmerHMMChan mode flatten
  file "*-dialign.fa"
  file "*-dialign.hmm"
  file "*-dialign.cons.fa"
  file "*-dialign.trimmed-cons.fa"
  file "*-dialign.trimmed.hmm"
  file "*-dialign.cons.vs_refmsacons"
  file "*-dialign.trimmed-cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run Dialign
  #export DIALIGN2_DIR=${dialignDir}/dialign2_dir
  #${dialignDir}/dialign2-2 -n -fa -fn ${simPrefix}-dialign ${referenceSeqFile} 
  #### Run Dialigntx
  ${dialignDir}/bin/dialign-tx -D ${dialignDir}/conf ${referenceSeqFile} ${simPrefix}-dialign.fa


  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-dialign.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-dialign.fa > ${simPrefix}-dialign.cons.fa
  ${projDir}/util/trimUnalignedEdges.pl ${simPrefix}-dialign.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-dialign.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-dialign.cons.fa L2.fa > gput100-train-dialign.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-dialign.trimmed-cons.fa L2.fa > gput100-train-dialign.trimmed-cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-dialign.cons.fa ${referenceMSACons} > ${simPrefix}-dialign.cons.vs_refmsacons 
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-dialign.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-dialign.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-dialign.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-dialign.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-dialign.hmm normal.stk > hmmbuild.log
  """
}

process runClustalW2 {
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForClustalW2

  when:
  runClustalW2

  output:
  tuple file("*clustalw2.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into clustalw2ToAMAChan
  tuple file("*clustalw2.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into clustalw2ToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into clustalw2ToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into clustalw2ToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) into clustalw2ToNhmmerHMMChan mode flatten
  file "*-clustalw2.fa"
  file "*-clustalw2.hmm"
  file "*-clustalw2.cons.fa"
  file "*-clustalw2.trimmed-cons.fa"
  file "*-clustalw2.trimmed.hmm"
  file "*-clustalw2.cons.vs_refmsacons"
  file "*-clustalw2.trimmed-cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run ClustalW2
  /usr/local/bin/clustalw2 -infile=${referenceSeqFile} -align -outfile=mangled.fa -output=FASTA
  cat mangled.fa | perl -ne '{ if ( /^>node-\\d+/ ) { s/_/:/; } print; }' > ${simPrefix}-clustalw2.fa 

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-clustalw2.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-clustalw2.fa > ${simPrefix}-clustalw2.cons.fa
  ${projDir}/util/trimUnalignedEdges.pl ${simPrefix}-clustalw2.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-clustalw2.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-clustalw2.cons.fa L2.fa > gput100-train-clustalw2.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-clustalw2.trimmed-cons.fa L2.fa > gput100-train-clustalw2.trimmed-cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-clustalw2.cons.fa ${referenceMSACons} > ${simPrefix}-clustalw2.cons.vs_refmsacons 
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-clustalw2.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-clustalw2.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-clustalw2.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-clustalw2.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-clustalw2.hmm normal.stk > hmmbuild.log
  """
}

process runClustalO {
  cpus 16
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForClustalO

  when:
  runClustalO

  output:
  tuple file("*clustalo.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into clustaloToAMAChan
  tuple file("*clustalo.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into clustaloToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into clustaloToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into clustaloToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) into clustaloToNhmmerHMMChan mode flatten
  file "*-clustalo.fa"
  file "*-clustalo.hmm"
  file "*-clustalo.cons.fa"
  file "*-clustalo.trimmed-cons.fa"
  file "*-clustalo.trimmed.hmm"
  file "*-clustalo.cons.vs_refmsacons"
  file "*-clustalo.trimmed-cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run ClustalOmega
  ${clustalOmegaDir}/clustalo --infile=${referenceSeqFile} --outfile=${simPrefix}-clustalo.fa --outfmt=FASTA --threads 16
  ## This doesn't seem to be a problem with clustalomega
  ##${clustalOmegaDir}/clustalo --infile=${referenceSeqFile} --outfile=mangled.fa --outfmt=FASTA
  ##cat mangled.fa | perl -ne '{ if ( /^>node-\\d+/ ) { s/_/:/; } print; }' > ${simPrefix}-clustalo.fa 

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-clustalo.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-clustalo.fa > ${simPrefix}-clustalo.cons.fa
  ${projDir}/util/trimUnalignedEdges.pl ${simPrefix}-clustalo.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-clustalo.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-clustalo.cons.fa L2.fa > gput100-train-clustalo.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-clustalo.trimmed-cons.fa L2.fa > gput100-train-clustalo.trimmed-cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-clustalo.cons.fa ${referenceMSACons} > ${simPrefix}-clustalo.cons.vs_refmsacons 
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-clustalo.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-clustalo.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-clustalo.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-clustalo.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-clustalo.hmm normal.stk > hmmbuild.log
  """
}

process runOpal {
  cpus 2
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForOpal

  when:
  runOpal

  output:
  tuple file("*opal.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into opalToAMAChan
  tuple file("*opal.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into opalToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into opalToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into opalToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) into opalToNhmmerHMMChan mode flatten
  file "*-opal.fa"
  file "*-opal.hmm"
  file "*-opal.cons.fa"
  file "*-opal.trimmed-cons.fa"
  file "*-opal.trimmed.hmm"
  file "*-opal.cons.vs_refmsacons"
  file "*-opal.trimmed-cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run Opal
  ${opalDir}/opal --mem 4G ${referenceSeqFile} > ${simPrefix}-opal.fa

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-opal.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-opal.fa > ${simPrefix}-opal.cons.fa
  ${projDir}/util/trimUnalignedEdges.pl ${simPrefix}-opal.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-opal.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-opal.cons.fa L2.fa > gput100-train-opal.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-opal.trimmed-cons.fa L2.fa > gput100-train-opal.trimmed-cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-opal.cons.fa ${referenceMSACons} > ${simPrefix}-opal.cons.vs_refmsacons 
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-opal.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-opal.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-opal.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-opal.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-opal.hmm normal.stk > hmmbuild.log
  """
}


process runMafft {
  cpus 16
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForMafft

  when:
  runMafft

  output:
  tuple file("*mafft.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into mafftToAMAChan
  tuple file("*mafft.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into mafftToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into mafftToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into mafftToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) into mafftToNhmmerHMMChan mode flatten
  file "*-mafft.fa"
  file "*-mafft.hmm"
  file "*-mafft.cons.fa"
  file "*-mafft.trimmed-cons.fa"
  file "*-mafft.trimmed.hmm"
  file "*-mafft.cons.vs_refmsacons"
  file "*-mafft.trimmed-cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run MAFFT
  # generate a filename like gput100-train-mafft.fa for output
  # E-INS-i :  --genafpair --maxiterate 1000
  # L-INS-i :  --localpair --maxiterate 1000
  # G-INS-i :  --globalpair --maxiterate 1000
  ${mafftDir}/mafft --thread 16 --localpair --maxiterate 1000 ${referenceSeqFile} | perl -ne '{ if ( /^>/ ) { print; } else { print uc(\$_); } }' > ${simPrefix}-mafft.fa 

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-mafft.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-mafft.fa > ${simPrefix}-mafft.cons.fa
  ${projDir}/util/trimUnalignedEdges.pl ${simPrefix}-mafft.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-mafft.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-mafft.cons.fa L2.fa > gput100-train-mafft.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-mafft.trimmed-cons.fa L2.fa > gput100-train-mafft.trimmed-cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-mafft.cons.fa ${referenceMSACons} > ${simPrefix}-mafft.cons.vs_refmsacons 
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-mafft.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-mafft.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-mafft.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-mafft.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-mafft.hmm normal.stk > hmmbuild.log
  """
}


process runMuscle {
  cpus 2
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForMuscle

  when:
  runMuscle

  output:
  tuple file("*muscle.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into muscleToAMAChan
  tuple file("*muscle.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into muscleToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into muscleToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into muscleToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) into muscleToNhmmerHMMChan mode flatten
  file "*-muscle.fa"
  file "*-muscle.hmm"
  file "*-muscle.cons.fa"
  file "*-muscle.trimmed.hmm"
  file "*-muscle.trimmed-cons.fa"
  file "*-muscle.cons.vs_refmsacons"
  file "*-muscle.trimmed-cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  # Run muscle and generate a filename like gput100-train-muscle.fa for output
  ${muscleDir}/muscle -in ${referenceSeqFile} -out ${simPrefix}-muscle.fa
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-muscle.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-muscle.fa > ${simPrefix}-muscle.cons.fa
  ${projDir}/util/trimUnalignedEdges.pl ${simPrefix}-muscle.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-muscle.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-muscle.cons.fa L2.fa > gput100-train-muscle.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-muscle.trimmed-cons.fa L2.fa > gput100-train-muscle.trimmed-cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-muscle.cons.fa ${referenceMSACons} > ${simPrefix}-muscle.cons.vs_refmsacons 
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-muscle.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-muscle.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-muscle.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-muscle.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-muscle.hmm normal.stk > hmmbuild.log
  """
}


process runRefiner {
  cpus 2
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(referenceMSACons), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from benchmarkFilesForRefiner

  when:
  runRefiner

  output:
  tuple file("*refiner-padded.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into refinerToAMAChan
  tuple file("*refiner-padded.fa"), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) into refinerToQScoreChan
  tuple file("*cons.fa"), file(testSeqFile) into refinerToCrossmatchChan mode flatten
  tuple file("*cons.fa"), file(testSeqFile) into refinerToNhmmerCONSChan mode flatten
  tuple file("*.hmm"), file(testSeqFile) optional true into refinerToNhmmerHMMChan mode flatten
  file "*refiner.log"
  file "*-refiner.stk"
  file "*-refiner.fa"
  file "*-refiner.hmm" optional true
  file "*-refiner.cons.fa"
  file "*-refiner.cons.vs_refmsacons"
  file "*-refiner-padded.hmm" optional true
  file "*-refiner-padded.cons.fa"
  file "*-refiner-padded.cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /^(.*-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  # Run refiner and generate a filename like gput100-train-refiner.fa for output
  ${repeatmodelerDir}/Refiner ${referenceSeqFile} >& ${simPrefix}-refiner.log
  ## eval Auto Run Blocker
  #${repeatmodelerDir}/util/Linup ${referenceSeqFile}.refiner.stk > alistart
  #${repeatmodelerDir}/util/AutoRunBlocker.pl -l alistart -w 7 -mc 4 -mr 2 > cons
  #${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  #${repeatmodelerDir}/util/AutoRunBlocker.pl -l alistart -w 15 -mc 4 -mr 2 > cons
  #${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  #${repeatmodelerDir}/util/AutoRunBlocker.pl -l alistart -w 24 -mc 4 -mr 2 > cons
  #${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  #${repeatmodelerDir}/util/AutoRunBlocker.pl -l alistart -w 5 -mc 4 -mr 2 > cons
  #${repeatmodelerDir}/util/alignAndCallConsensus.pl -re -c cons -e ${referenceSeqFile}
  #${repeatmodelerDir}/util/Linup -stockholm -name ${referenceMSAFile.baseName} rep.out > ${simPrefix}-refiner.stk
  #mv cons ${simPrefix}-refiner.cons.fa
  ## End Auto Run Blocker
  # Generates *.refiner.stk and *.refiner_cons ... rename to final files
  mv ${referenceSeqFile}.refiner.stk ${simPrefix}-refiner.stk
  mv ${referenceSeqFile}.refiner_cons ${simPrefix}-refiner.cons.fa
  ${repeatmodelerDir}/util/Linup -msa ${simPrefix}-refiner.stk > ${simPrefix}-refiner.fa
  ${projDir}/util/stkToQscoreMSA.pl ${referenceSeqFile} ${simPrefix}-refiner.stk > ${simPrefix}-refiner-padded.fa
  ${projDir}/util/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-refiner-padded.fa

  #### Generate Consensus Model
  # Generate a consensus from the padded alignment
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-refiner-padded.fa > ${simPrefix}-refiner-padded.cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-refiner-padded.cons.fa L2.fa > gput100-train-refiner-padded.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-refiner.cons.fa L2.fa > gput100-train-refiner.cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-refiner-padded.cons.fa ${referenceMSACons} > ${simPrefix}-refiner-padded.cons.vs_refmsacons
  ${projDir}/util/compareConsensiNeedle.pl ${simPrefix}-refiner.cons.fa ${referenceMSACons} > ${simPrefix}-refiner.cons.vs_refmsacons
 
  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-refiner-padded.fa > padded.stk
  if ! ${hmmerDir}/hmmbuild ${simPrefix}-refiner-padded.hmm padded.stk > padded-hmmbuild.log; then
    echo "No padded hmm"
    if [ -f ${simPrefix}-refiner-padded.hmm ]; then
      rm ${simPrefix}-refiner-padded.hmm
    fi
  fi
  if ! ${hmmerDir}/hmmbuild ${simPrefix}-refiner.hmm ${simPrefix}-refiner.stk > hmmbuild.log; then
    echo "No hmm"
    if [ -f ${simPrefix}-refiner.hmm ]; then
      rm ${simPrefix}-refiner.hmm
    fi
  fi
  """
}

process evalAMA {
  cpus 2
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  when:
  runDartScore

  input:
  set file(predictedMSAFile), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from refinerToAMAChan.mix(muscleToAMAChan,mafftToAMAChan,clustalw2ToAMAChan,dialignToAMAChan,kalignToAMAChan,fsaToAMAChan, opalToAMAChan, clustaloToAMAChan)

  output:
  file "*.ama_score"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceSeqFile.toRealPath().getName(referenceSeqFile.toRealPath().getNameCount() - 2)
  """
  ${projDir}/util/fastaMSAToSTK.pl ${predictedMSAFile} > ${predictedMSAFile.baseName}.stk
  ${projDir}/util/fastaMSAToSTK.pl ${referenceMSAFile} > reference.stk
  ${dartDir}/cmpalign reference.stk ${predictedMSAFile.baseName}.stk >& ${predictedMSAFile.baseName}.ama_score
  """
}

process evalQScore {
  cpus 2
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(predictedMSAFile), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from refinerToQScoreChan.mix(muscleToQScoreChan,mafftToQScoreChan,clustalw2ToQScoreChan,dialignToQScoreChan,kalignToQScoreChan,fsaToQScoreChan,opalToQScoreChan,clustaloToQScoreChan)

  when:
  runQScore

  output:
  file "*.qscore_score"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceSeqFile.toRealPath().getName(referenceSeqFile.toRealPath().getNameCount() - 2)
  """
  ${qscoreDir}/qscore  -test ${predictedMSAFile} -ref ${referenceMSAFile} >& ${predictedMSAFile.baseName}.qscore_score
  """
}


process runCrossmatch {
  cpus 2
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  when:
  runCMConsEval

  input:
  set file(consFile), file(testSeqFile) from refinerToCrossmatchChan.mix(muscleToCrossmatchChan,mafftToCrossmatchChan,clustalw2ToCrossmatchChan,dialignToCrossmatchChan,kalignToCrossmatchChan,fsaToCrossmatchChan,opalToCrossmatchChan,clustaloToCrossmatchChan)

  output:
  file "*.cm_score"
  file "*.cm"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = testSeqFile.toRealPath().getName(testSeqFile.toRealPath().getNameCount() - 2)
  """
  ${phrapDir}/cross_match -matrix ${projDir}/matrices/14p35g.matrix -masklevel 80 -gap_init -30 -ins_gap_ext -6 -del_gap_ext -5 -minmatch 7 -minscore 180 ${testSeqFile} ${consFile}  > ${consFile.baseName}.cm
  cat ${consFile.baseName}.cm | perl -ne '{if ( /^\\s*(\\d+)\\s+\\d+\\.\\d+/ ){ \$sum+=\$1; } print \$sum if ( eof );}' > ${consFile.baseName}.cm_score
  """
}


process runNhmmerHMM {
  cpus 8
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  when:
  runNhmmerHMMEval

  input:
  set file(hmmFile), file(testSeqFile) from refinerToNhmmerHMMChan.mix(muscleToNhmmerHMMChan,mafftToNhmmerHMMChan,clustalw2ToNhmmerHMMChan,dialignToNhmmerHMMChan,kalignToNhmmerHMMChan,fsaToNhmmerHMMChan,opalToNhmmerHMMChan,clustaloToNhmmerHMMChan)

  output:
  file "*nhmmer"
  file "*nhmmer_score"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = testSeqFile.toRealPath().getName(testSeqFile.toRealPath().getNameCount() - 2)
  """
  ${hmmerDir}/nhmmer --cpu 8 --noali --dfamtblout ${hmmFile.baseName}.nhmmer ${hmmFile} ${testSeqFile} > /dev/null
  cat ${hmmFile.baseName}.nhmmer | tr -s ' ' | cut -d ' ' -f 4 | awk '{s+=\$1} END {print s}' > ${hmmFile.baseName}.nhmmer_score
  """
}

process runNhmmerCONS {
  cpus 8
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  when:
  runNhmmerConsEval

  input:
  set file(consFile), file(testSeqFile) from refinerToNhmmerCONSChan.mix(muscleToNhmmerCONSChan,mafftToNhmmerCONSChan,clustalw2ToNhmmerCONSChan,dialignToNhmmerCONSChan,kalignToNhmmerCONSChan,fsaToNhmmerCONSChan,opalToNhmmerCONSChan,clustaloToNhmmerCONSChan)

  output:
  file "*nhmmer"
  file "*nhmmer_score"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = testSeqFile.toRealPath().getName(testSeqFile.toRealPath().getNameCount() - 2)
  """
  ${hmmerDir}/nhmmer --cpu 8 --dna --noali --dfamtblout ${consFile.baseName}.nhmmer ${consFile} ${testSeqFile} > /dev/null
  cat ${consFile.baseName}.nhmmer | tr -s ' ' | cut -d ' ' -f 4 | awk '{s+=\$1} END {print s}' > ${consFile.baseName}.nhmmer_score
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


// ATTIC
//
//process evalFastSP {
//  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }
//
//  input:
//  set file(predictedMSAFile), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from predictedRefinerChan3.mix(predictedMuscleChan3,predictedMafftChan3)
//
//  when:
//  runFastSP
//
//  output:
//  file "*.fastsp_score"
//
//  script:
//  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
//  repDir = referenceSeqFile.toRealPath().getName(referenceSeqFile.toRealPath().getNameCount() - 2)
//  """
//  java -jar ${fastSPDir}/FastSP.jar -r ${referenceMSAFile} -e ${predictedMSAFile} >& ${predictedMSAFile.baseName}.fastsp_score
//  """
//}
//


