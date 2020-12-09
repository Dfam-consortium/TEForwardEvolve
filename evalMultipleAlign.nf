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
     --dialign            : Run Dialign [optional]
     --kalign             : Run Kalign [optional]
     --clustalw2          : Run ClustalW2 [optional]
     --mafft              : Run Mafft [optiona]
     --cluster            : Either "local", "quanah", "hrothgar" or "griz" 
                            default="local"

 Others to evaluate:

  Prank: 
    /usr/local/prank/bin/prank -d=ff.seqs
    Took 85 minutes for one alignment!
    Changes ":" to "_"

  T-Coffee: 
     /nfs/public/ro/es/appbin/linux-x86_64/T-COFFEE_installer_Version_13.41.0.28bdc39_linux_x64/bin/t_coffee -in tcoffee-E20201208-235308-0480-65402959-p1m.sequence -case=upper -n_core=8 -output=clustalw,msf,phylip,score_html,fasta -outorder=aligned -type=dna; echo ' '
    /usr/local/src/T-COFFEE_installer_Version_13.45.0.4846264_linux_x64/bin/t_coffee -in ff.seqs -case=upper -n_core=8 -output=fasta -outorder=aligned -type=dna
    # uses 75G of RAM before it was killed  ( run at EBI took 20 minutes...at tcoffee.crg.cat it took 2 hrs 7 min )
    # Changes ":" to "_"
  
  Clustal Omega:

    
Test=gput2750-train-clustalw2.fa;Ref=gput2750-train-refmsa.fa;Q=0.0307;TC=0.00585
Test=gput2750-train-mafft.fa;Ref=gput2750-train-refmsa.fa;Q=0.695;TC=0.0406
Test=gput2750-train-muscle.fa;Ref=gput2750-train-refmsa.fa;Q=0.346;TC=0.0143
Test=gput2750-train-refiner-padded.fa;Ref=gput2750-train-refmsa.fa;Q=0.669;TC=0.0179

Test=dialign-msa.fa;Ref=rep-1/gput2750-train-refmsa.fa;Q=0.515;TC=0.0227
Test=prank-msa.fa;Ref=rep-1/gput2750-train-refmsa.fa;Q=0.362;TC=0.0117
Test=tcoffee-msa.fa;Ref=rep-1/gput2750-train-refmsa.fa;Q=0.3;TC=0.00622
Test=kalign-msa.fa;Ref=rep-1/gput2750-train-refmsa.fa;Q=0.23;TC=0.00512


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
params.dialign = false
params.kalign = false
params.benchmarkDir = "${workflow.projectDir}/sample"
params.outputDir = "results"
params.seed = "${workflow.projectDir}/sample/L2.fa"

runMuscle = params.muscle
runRefiner = params.refiner
runMafft = params.mafft
runClustalW2 = params.clustalw2
runDialign = params.dialign
runKalign = params.kalign
runFastSP = false
runQScore = true
outputDir = params.outputDir

// Default software dependencies ( see localizations in cluster sections )
qscoreDir = "/home/rhubley/projects/DNAMultipleAlignment/qscore"
phrapDir = "/usr/local/phrap"
hmmerDir = "/usr/local/hmmer/bin"
mafftDir = "/usr/local/mafft/bin"
dialignDir = "/usr/local/dialign-2.2.1"
kalignDir = "/usr/local/kalign2"
clustalW2Dir = "/usr/local/bin"
dartDir = "/home/rhubley/projects/DNAMultipleAlignment/dart/bin"
muscleDir = "/usr/local/bin"
repeatmodelerDir = "/home/rhubley/projects/RepeatModeler"
fastSPDir = "/home/rhubley/projects/DNAMultipleAlignment/FastSP"

//FOR DEBUGGING...limit the files run
//Channel.fromFilePairs( params.benchmarkDir + "/rep-1/gput100-{train-seqs,train-refmsa,test-seqs}.fa", size: 3, flat:true )

//
// E.g "gput100-", "gput100-train-seqs.fa", "gput100-train-refmsa.fa", "gput100-test-seqs.fa"
Channel.fromFilePairs( params.benchmarkDir + "/*/gput*-{train-seqs,train-refmsa,test-seqs}.fa", size: 3, flat:true )
       .into { benchmarkFilesForComp }

//       .into { benchmarkFilesForComp; benchmarkFilesForMafft; benchmarkFilesForRefiner; benchmarkFilesForMuscle; benchmarkFilesForClustalW2 }

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

log.info "evalMultipleAlign.nf : Multiple Alignment Evaluation ver 0.2"
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
log.info "Cluster             : " + params.cluster
log.info "Output DIR          : " + outputDir
log.info "Queue/Partititon    : " + thisQueue






process processRefMSA {
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
  file "*cons.fa"
  file "*cons.vs_self"
  file "*cons.vs_seed"
  file "*avgKDiv"
  

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /(gput\d+-train)/)[0][0]
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
  ${workflow.projectDir}/compareConsensiNeedle.pl ${referenceMSAFile.baseName}.cons.fa ${seedFile} > ${referenceMSAFile.baseName}.cons.vs_seed
  ${workflow.projectDir}/compareConsensiNeedle.pl ${referenceMSAFile.baseName}.cons.fa ${referenceMSAFile.baseName}.cons.fa > ${referenceMSAFile.baseName}.cons.vs_self
  ###########
  """
}


//  KAlign: https://msa.sbc.su.se/cgi-bin/msa.cgi
//    parameters suggested by website for DNA
process runKalign {
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
  simPrefix = (referenceMSAFile.name =~ /(gput\d+-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run Kalign
  ${kalignDir}/kalign -gpo 80 -gpe 3 -tgpe 3 -bonus 0 ${referenceSeqFile} > ${simPrefix}-kalign.fa

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${workflow.projectDir}/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-kalign.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-kalign.fa > ${simPrefix}-kalign.cons.fa
  ${workflow.projectDir}/trimUnalignedEdges.pl ${simPrefix}-kalign.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-kalign.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-kalign.cons.fa L2.fa > gput100-train-kalign.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-kalign.trimmed-cons.fa L2.fa > gput100-train-kalign.trimmed-cons.vs_refmsacons
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-kalign.cons.fa ${referenceMSACons} > ${simPrefix}-kalign.cons.vs_refmsacons 
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-kalign.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-kalign.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-kalign.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-kalign.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-kalign.hmm normal.stk > hmmbuild.log
  """
}

//  DIALIGN: http://dialign.gobics.de/
process runDialign {
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
  simPrefix = (referenceMSAFile.name =~ /(gput\d+-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run Dialign
  export DIALIGN2_DIR=${dialignDir}/dialign2_dir
  ${dialignDir}/dialign2-2 -n -fa -fn ${simPrefix}-dialign ${referenceSeqFile} 

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${workflow.projectDir}/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-dialign.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-dialign.fa > ${simPrefix}-dialign.cons.fa
  ${workflow.projectDir}/trimUnalignedEdges.pl ${simPrefix}-dialign.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-dialign.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-dialign.cons.fa L2.fa > gput100-train-dialign.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-dialign.trimmed-cons.fa L2.fa > gput100-train-dialign.trimmed-cons.vs_refmsacons
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-dialign.cons.fa ${referenceMSACons} > ${simPrefix}-dialign.cons.vs_refmsacons 
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-dialign.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-dialign.trimmed-cons.vs_refmsacons

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
  simPrefix = (referenceMSAFile.name =~ /(gput\d+-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run ClustalW2
  /usr/local/bin/clustalw2 -infile=${referenceSeqFile} -align -outfile=mangled.fa -output=FASTA
  cat mangled.fa | perl -ne '{ if ( /^>node-\\d+/ ) { s/_/:/; } print; }' > ${simPrefix}-clustalw2.fa 

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${workflow.projectDir}/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-clustalw2.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-clustalw2.fa > ${simPrefix}-clustalw2.cons.fa
  ${workflow.projectDir}/trimUnalignedEdges.pl ${simPrefix}-clustalw2.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-clustalw2.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-clustalw2.cons.fa L2.fa > gput100-train-clustalw2.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-clustalw2.trimmed-cons.fa L2.fa > gput100-train-clustalw2.trimmed-cons.vs_refmsacons
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-clustalw2.cons.fa ${referenceMSACons} > ${simPrefix}-clustalw2.cons.vs_refmsacons 
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-clustalw2.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-clustalw2.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-clustalw2.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-clustalw2.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-clustalw2.hmm normal.stk > hmmbuild.log
  """
}

 

process runMafft {
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
  simPrefix = (referenceMSAFile.name =~ /(gput\d+-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  #### Run MAFFT
  # generate a filename like gput100-train-mafft.fa for output
  # E-INS-i :  --genafpair --maxiterate 1000
  # L-INS-i :  --localpair --maxiterate 1000
  # G-INS-i :  --globalpair --maxiterate 1000
  ${mafftDir}/mafft --localpair --maxiterate 1000 ${referenceSeqFile} | perl -ne '{ if ( /^>/ ) { print; } else { print uc(\$_); } }' > ${simPrefix}-mafft.fa 

  #### Sanity Check MSA
  # Since this is a full sequence MSA the validate the sequences against the reference MSA as a sanity check
  ${workflow.projectDir}/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-mafft.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-mafft.fa > ${simPrefix}-mafft.cons.fa
  ${workflow.projectDir}/trimUnalignedEdges.pl ${simPrefix}-mafft.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-mafft.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-mafft.cons.fa L2.fa > gput100-train-mafft.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-mafft.trimmed-cons.fa L2.fa > gput100-train-mafft.trimmed-cons.vs_refmsacons
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-mafft.cons.fa ${referenceMSACons} > ${simPrefix}-mafft.cons.vs_refmsacons 
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-mafft.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-mafft.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-mafft.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-mafft.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-mafft.hmm normal.stk > hmmbuild.log
  """
}


process runMuscle {
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
  simPrefix = (referenceMSAFile.name =~ /(gput\d+-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  # Run muscle and generate a filename like gput100-train-muscle.fa for output
  ${muscleDir}/muscle -in ${referenceSeqFile} -out ${simPrefix}-muscle.fa
  ${workflow.projectDir}/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-muscle.fa

  #### Generate Consensus Model
  # Generate two consensus files from the multiple alignment.  One which includes single sequence alignment edges, and a trimmed version
  # which trims back until there is at least one column containing two or more sequences.
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-muscle.fa > ${simPrefix}-muscle.cons.fa
  ${workflow.projectDir}/trimUnalignedEdges.pl ${simPrefix}-muscle.fa > trimmed-msa.fa
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} trimmed-msa.fa > ${simPrefix}-muscle.trimmed-cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-muscle.cons.fa L2.fa > gput100-train-muscle.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-muscle.trimmed-cons.fa L2.fa > gput100-train-muscle.trimmed-cons.vs_refmsacons
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-muscle.cons.fa ${referenceMSACons} > ${simPrefix}-muscle.cons.vs_refmsacons 
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-muscle.trimmed-cons.fa ${referenceMSACons} > ${simPrefix}-muscle.trimmed-cons.vs_refmsacons

  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm trimmed-msa.fa > trimmed.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-muscle.trimmed.hmm trimmed.stk > trimmed-hmmbuild.log
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-muscle.fa > normal.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-muscle.hmm normal.stk > hmmbuild.log
  """
}


process runRefiner {
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
  tuple file("*.hmm"), file(testSeqFile) into refinerToNhmmerHMMChan mode flatten
  file "*refiner.log"
  file "*-refiner.stk"
  file "*-refiner.fa"
  file "*-refiner.hmm"
  file "*-refiner.cons.fa"
  file "*-refiner.cons.vs_refmsacons"
  file "*-refiner-padded.hmm"
  file "*-refiner-padded.cons.fa"
  file "*-refiner-padded.cons.vs_refmsacons"

  script:
  // Identify the prefix "gput100-train" from the filename for use in the script
  simPrefix = (referenceMSAFile.name =~ /(gput\d+-train)/)[0][0]
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceMSAFile.toRealPath().getName(referenceMSAFile.toRealPath().getNameCount() - 2)
  """
  # Run refiner and generate a filename like gput100-train-refiner.fa for output
  ${repeatmodelerDir}/Refiner ${referenceSeqFile} >& ${simPrefix}-refiner.log
  ## eval Auto Run Blocker
  ${repeatmodelerDir}/util/Linup ${referenceSeqFile}.refiner.stk > alistart
  ${repeatmodelerDir}/util/protocol/AutoRunBlocker.pl -l alistart -w 7 -mc 4 -mr 2 > cons
  ${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  ${repeatmodelerDir}/util/protocol/AutoRunBlocker.pl -l alistart -w 15 -mc 4 -mr 2 > cons
  ${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  ${repeatmodelerDir}/util/protocol/AutoRunBlocker.pl -l alistart -w 24 -mc 4 -mr 2 > cons
  ${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  ${repeatmodelerDir}/util/protocol/AutoRunBlocker.pl -l alistart -w 5 -mc 4 -mr 2 > cons
  ${repeatmodelerDir}/util/alignAndCallConsensus.pl -re -c cons -e ${referenceSeqFile}
  ${repeatmodelerDir}/util/Linup -stockholm -name ${referenceMSAFile.baseName} rep.out > ${simPrefix}-refiner.stk
  mv cons ${simPrefix}-refiner.cons.fa
  ## End Auto Run Blocker
  # Generates *.refiner.stk and *.refiner_cons ... rename to final files
  #mv ${referenceSeqFile}.refiner.stk ${simPrefix}-refiner.stk
  #mv ${referenceSeqFile}.refiner_cons ${simPrefix}-refiner.cons.fa
  ${repeatmodelerDir}/util/Linup -msa ${simPrefix}-refiner.stk > ${simPrefix}-refiner.fa
  ${workflow.projectDir}/stkToQscoreMSA.pl ${referenceSeqFile} ${simPrefix}-refiner.stk > ${simPrefix}-refiner-padded.fa
  ${workflow.projectDir}/validateEstimatedMSA.pl ${referenceMSAFile} ${simPrefix}-refiner-padded.fa

  #### Generate Consensus Model
  # Generate a consensus from the padded alignment
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceMSAFile.baseName} ${simPrefix}-refiner-padded.fa > ${simPrefix}-refiner-padded.cons.fa

  #### Consensus VS Consensus
  # Compare consensi to evalute how well each can rebuild the ref msa consensus
  #   E.g. compareConsensiNeedle.pl gput100-train-refiner-padded.cons.fa L2.fa > gput100-train-refiner-padded.cons.vs_refmsacons
  #        compareConsensiNeedle.pl gput100-train-refiner.cons.fa L2.fa > gput100-train-refiner.cons.vs_refmsacons
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-refiner-padded.cons.fa ${referenceMSACons} > ${simPrefix}-refiner-padded.cons.vs_refmsacons
  ${workflow.projectDir}/compareConsensiNeedle.pl ${simPrefix}-refiner.cons.fa ${referenceMSACons} > ${simPrefix}-refiner.cons.vs_refmsacons
 
  #### Generate HMM Model
  ${repeatmodelerDir}/util/Linup -noTemplate -name predicted -stockholm ${simPrefix}-refiner-padded.fa > padded.stk
  ${hmmerDir}/hmmbuild ${simPrefix}-refiner-padded.hmm padded.stk > padded-hmmbuild.log
  ${hmmerDir}/hmmbuild ${simPrefix}-refiner.hmm ${simPrefix}-refiner.stk > hmmbuild.log
  """
}

process evalAMA {
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(predictedMSAFile), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from refinerToAMAChan.mix(muscleToAMAChan,mafftToAMAChan,clustalw2ToAMAChan,dialignToAMAChan,kalignToAMAChan)

  output:
  file "*.ama_score"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = referenceSeqFile.toRealPath().getName(referenceSeqFile.toRealPath().getNameCount() - 2)
  """
  ${workflow.projectDir}/fastaMSAToSTK.pl ${predictedMSAFile} > ${predictedMSAFile.baseName}.stk
  ${workflow.projectDir}/fastaMSAToSTK.pl ${referenceMSAFile} > reference.stk
  ${dartDir}/cmpalign reference.stk ${predictedMSAFile.baseName}.stk >& ${predictedMSAFile.baseName}.ama_score
  """
}

process evalQScore {
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(predictedMSAFile), file(testSeqFile), file(referenceMSAFile), file(referenceSeqFile) from refinerToQScoreChan.mix(muscleToQScoreChan,mafftToQScoreChan,clustalw2ToQScoreChan,dialignToQScoreChan,kalignToQScoreChan)

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
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(consFile), file(testSeqFile) from refinerToCrossmatchChan.mix(muscleToCrossmatchChan,mafftToCrossmatchChan,clustalw2ToCrossmatchChan,dialignToCrossmatchChan,kalignToCrossmatchChan)

  output:
  file "*.cm_score"
  file "*.cm"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = testSeqFile.toRealPath().getName(testSeqFile.toRealPath().getNameCount() - 2)
  """
  ${phrapDir}/cross_match -matrix ${workflow.projectDir}/matrices/14p35g.matrix -masklevel 80 -gap_init -30 -ins_gap_ext -6 -del_gap_ext -5 -minmatch 7 -minscore 180 ${testSeqFile} ${consFile}  > ${consFile.baseName}.cm
  cat ${consFile.baseName}.cm | perl -ne '{if ( /^\\s*(\\d+)\\s+\\d+\\.\\d+/ ){ \$sum+=\$1; } print \$sum if ( eof );}' > ${consFile.baseName}.cm_score
  """
}


process runNhmmerHMM {
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(hmmFile), file(testSeqFile) from refinerToNhmmerHMMChan.mix(muscleToNhmmerHMMChan,mafftToNhmmerHMMChan,clustalw2ToNhmmerHMMChan,dialignToNhmmerHMMChan,kalignToNhmmerHMMChan)

  output:
  file "*nhmmer"
  file "*nhmmer_score"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = testSeqFile.toRealPath().getName(testSeqFile.toRealPath().getNameCount() - 2)
  """
  ${hmmerDir}/nhmmer --noali --dfamtblout ${hmmFile.baseName}.nhmmer ${hmmFile} ${testSeqFile} > /dev/null
  cat ${hmmFile.baseName}.nhmmer | tr -s ' ' | cut -d ' ' -f 4 | awk '{s+=\$1} END {print s}' > ${hmmFile.baseName}.nhmmer_score
  """
}

process runNhmmerCONS {
  publishDir "${outputDir}", mode: 'copy', saveAs: { filename -> "$repDir/$filename" }

  input:
  set file(consFile), file(testSeqFile) from refinerToNhmmerCONSChan.mix(muscleToNhmmerCONSChan,mafftToNhmmerCONSChan,clustalw2ToNhmmerCONSChan,dialignToNhmmerCONSChan,kalignToNhmmerCONSChan)

  output:
  file "*nhmmer"
  file "*nhmmer_score"

  script:
  // Identify the "rep-#" directory from the path to the referenceMSAFile for output
  repDir = testSeqFile.toRealPath().getName(testSeqFile.toRealPath().getNameCount() - 2)
  """
  ${hmmerDir}/nhmmer --noali --dfamtblout ${consFile.baseName}.nhmmer ${consFile} ${testSeqFile} > /dev/null
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


