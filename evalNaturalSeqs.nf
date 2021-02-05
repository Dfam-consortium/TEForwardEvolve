#!/usr/bin/env nextflow
/*
vim: syntax=groovy

evalNaturalSeqs.nf : Evaluate multiple alignment programs on natural TE family

 Parameters:

     --outputDir <dir>    : Directory to store the results
     --familyDir <dir>    : Directory where a benchmark family is found
     --protein <file>     : File containing the hand-curated protein for the family
     --muscle             : Run Muscle [optional]
     --refiner            : Run Refiner [optional]
     --dialign            : Run DialignTX [optional]
     --kalign             : Run Kalign [optional]
     --fsa                : Run FSA [optional]
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

  FSA:

Robert Hubley, 12/2020
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
params.familyDir = "undefined"
params.outputDir = "results"
params.protein = "undefined"

runMuscle = params.muscle
runRefiner = params.refiner
runMafft = params.mafft
runClustalW2 = params.clustalw2
runDialign = params.dialign
runKalign = params.kalign
runFSA = params.fsa
outputDir = params.outputDir

// Default software dependencies ( see localizations in cluster sections )
blastDir = "/usr/local/rmblast/bin"
hmmerDir = "/usr/local/hmmer/bin"
mafftDir = "/usr/local/mafft/bin"
//dialignDir = "/usr/local/dialign-2.2.1"
dialignDir = "/usr/local/dialign-tx-1.0.2"
kalignDir = "/usr/local/kalign2"
clustalW2Dir = "/usr/local/bin"
fsaDir = "/usr/local/fsa/bin"
muscleDir = "/usr/local/bin"
exonerateDir = "/usr/local/exonerate-2.2.0-x86_64/bin"
repeatmodelerDir = "/home/rhubley/projects/RepeatModeler"

// Generate a series of samples from the orffrags file increasing in size:
//   100, 150, 200, 250, 300, 350, 400, 450, 500, 1000 
seqFile = file( params.familyDir + "/orffrags.fa")
//
proteinFile = file(params.protein)

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

log.info "evalNaturalSeqs.nf : Multiple Alignment Evaluation ver 0.1"
log.info "============================================================"
log.info "working directory   : " + workflow.workDir
log.info "RepeatModeler DIR   : " + repeatmodelerDir
log.info "Family DIR          : " + params.familyDir
log.info "Protein File        : " + params.protein
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
log.info "Cluster             : " + params.cluster
log.info "Output DIR          : " + outputDir
log.info "Queue/Partititon    : " + thisQueue



process generateSamples {
  publishDir "${outputDir}", mode: 'copy' 

  input:
  file seqFile
  
  output:
  file "*sample.fa" into benchmarkFilesForClustalW2 
  file "*sample.fa" into benchmarkFilesForMafft 
  file "*sample.fa" into benchmarkFilesForMuscle 
  file "*sample.fa" into benchmarkFilesForRefiner 
  file "*sample.fa" into benchmarkFilesForKalign 
  file "*sample.fa" into benchmarkFilesForDialign 
  file "*sample.fa" into benchmarkFilesForFSA 
  file "" into benchmarkFilesForBaseline
  file "*sample.fa"

  script:
  """
  # A hack to get the same sample sets over multiple runs of this Nextflow script:
  ${workflow.projectDir}/sampleFromFA.pl "100,150,200,250,300,350,400,450,500,1000" ${seqFile} 1608835988
  """
}


process runFSA {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file proteinFile
  file referenceSeqFile from benchmarkFilesForFSA.flatten()

  when:
  runFSA

  output:
  file "*-fsa.fa"
  file "*-fsa.cons.fa"
  file "*-fsa.blastx"
  file "*-fsa.exonerate"
  file "*-fsa.time"

  script:
  methodPrefix = "fsa"
  sampleSize = (referenceSeqFile.name =~ /^.*\-(\d+)sample.fa/)[0][1]
  //log.info("sampleSize = " + sampleSize)
  """
  if [ ${sampleSize} -lt 350 ] 
  then
    #### Run FSA
    SECONDS=0
    ${fsaDir}/fsa ${referenceSeqFile} > ${referenceSeqFile.baseName}-${methodPrefix}.fa
    echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time
  
    # Gen cons and eval
    ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  else 
    touch ${referenceSeqFile.baseName}-${methodPrefix}.fa
    touch ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
    echo "0" > ${referenceSeqFile.baseName}-${methodPrefix}.time
    touch ${referenceSeqFile.baseName}-${methodPrefix}.blastx
    touch ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  fi
  """
}

//  KAlign: https://msa.sbc.su.se/cgi-bin/msa.cgi
//    parameters suggested by website for DNA
process runKalign {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file proteinFile
  file referenceSeqFile from benchmarkFilesForKalign.flatten()

  when:
  runKalign

  output:
  file "*-kalign.fa"
  file "*-kalign.cons.fa"
  file "*-kalign.blastx"
  file "*-kalign.exonerate"
  file "*-kalign.time"

  script:
  methodPrefix = "kalign"
  """
  #### Run Kalign
  SECONDS=0
  ${kalignDir}/kalign -gpo 80 -gpe 3 -tgpe 3 -bonus 0 ${referenceSeqFile} > ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}

//  DIALIGN: http://dialign.gobics.de/
// Switched to dialigntx because it does appear to perform better
process runDialign {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file proteinFile
  file referenceSeqFile from benchmarkFilesForDialign.flatten()

  when:
  runDialign

  output:
  file "*-dialign.fa"
  file "*-dialign.cons.fa"
  file "*-dialign.blastx"
  file "*-dialign.exonerate"
  file "*-dialign.time"

  script:
  methodPrefix = "dialign"
  """
  #### Run Dialign
  #export DIALIGN2_DIR=${dialignDir}/dialign2_dir
  #${dialignDir}/dialign2-2 -n -fa -fn ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile} 
  #### Run Dialigntx
  SECONDS=0
  ${dialignDir}/dialign-tx -D ${dialignDir}/conf ${referenceSeqFile} ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time


  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}

process runClustalW2 {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file proteinFile
  file referenceSeqFile from benchmarkFilesForClustalW2.flatten()

  when:
  runClustalW2

  output:
  file "*-clustalw2.fa"
  file "*-clustalw2.cons.fa"
  file "*-clustalw2.blastx"
  file "*-clustalw2.exonerate"
  file "*-clustalw2.time"

  script:
  methodPrefix = "clustalw2"
  """
  #### Run ClustalW2
  SECONDS=0
  /usr/local/bin/clustalw2 -infile=${referenceSeqFile} -align -outfile=mangled.fa -output=FASTA
  cat mangled.fa | perl -ne '{ if ( /^>node-\\d+/ ) { s/_/:/; } print; }' > ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}

 

process runMafft {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file proteinFile
  file referenceSeqFile from benchmarkFilesForMafft.flatten()

  when:
  runMafft

  output:
  file "*-mafft.fa"
  file "*-mafft.cons.fa"
  file "*-mafft.blastx"
  file "*-mafft.exonerate"
  file "*-mafft.time"

  script:
  methodPrefix = "mafft"
  """
  #### Run MAFFT
  # generate a filename like gput100-train-mafft.fa for output
  # E-INS-i :  --genafpair --maxiterate 1000
  # L-INS-i :  --localpair --maxiterate 1000
  # G-INS-i :  --globalpair --maxiterate 1000
  SECONDS=0
  ${mafftDir}/mafft --localpair --maxiterate 1000 ${referenceSeqFile} | perl -ne '{ if ( /^>/ ) { print; } else { print uc(\$_); } }' > ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}


process runMuscle {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file proteinFile
  file referenceSeqFile from benchmarkFilesForMuscle.flatten()

  when:
  runMuscle

  output:
  file "*-muscle.fa"
  file "*-muscle.cons.fa"
  file "*-muscle.blastx"
  file "*-muscle.exonerate"
  file "*-muscle.time"

  script:
  methodPrefix = "muscle"
  """
  # Run muscle and generate a filename like gput100-train-muscle.fa for output
  SECONDS=0
  ${muscleDir}/muscle -in ${referenceSeqFile} -out ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}


process runRefiner {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file proteinFile
  file referenceSeqFile from benchmarkFilesForRefiner.flatten()

  when:
  runRefiner

  output:
  file "*-refiner.fa"
  file "*-refiner.cons.fa"
  file "*-refiner.blastx"
  file "*-refiner.exonerate"
  file "*-refiner.time"

  script:
  methodPrefix = "refiner"
  """
  SECONDS=0
  # Run refiner and generate a filename like gput100-train-refiner.fa for output
  ${repeatmodelerDir}/Refiner ${referenceSeqFile} >& ${referenceSeqFile.baseName}-${methodPrefix}.log
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
  ${repeatmodelerDir}/util/Linup -msa rep.out > ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${proteinFile} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
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
