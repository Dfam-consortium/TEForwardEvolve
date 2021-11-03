#!/usr/bin/env nextflow
/*
vim: syntax=groovy

evalNaturalSeqs.nf : Evaluate multiple alignment programs on natural TE family

 Parameters:

     --outputDir <dir>        : Directory to store the results
     --familyDir <dir>        : Directory where a benchmark family is found
     --curated_protein <file> : File containing the hand-curated protein for the family
     --db_proteins <file>     : File containing the related db proteins for the family
     --muscle                 : Run Muscle [optional]
     --refiner                : Run Refiner [optional]
     --dialign                : Run DialignTX [optional]
     --tcoffee                : Run TCoffee [optional]
     --probcons               : Run Probcons [optional]
     --kalign                 : Run Kalign [optional]
     --fsa                    : Run FSA [optional]
     --clustalo               : Run ClustalOmega [optional]
     --mafft                  : Run Mafft [optiona]
     --cluster                : Either "local", "quanah", "hrothgar" or "griz" 
                            default="local"

Robert Hubley, 12/2020
*/


// Defaults
params.cluster = "local"
params.muscle = false
params.refiner = false
params.clustalo = false
params.mafft = false
params.fsa = false
params.dialign = false
params.tcoffee = false
params.probcons = false
params.kalign = false
params.familyDir = "undefined"
params.outputDir = "results"
params.curated_protein = "undefined"
params.db_proteins = "undefined"

runMuscle = params.muscle
runRefiner = params.refiner
runMafft = params.mafft
runClustalO = params.clustalo
runDialign = params.dialign
runKalign = params.kalign
runTCoffee = params.tcoffee
runProbcons = params.probcons
runFSA = params.fsa
outputDir = params.outputDir

// Default software dependencies ( see localizations in cluster sections )
blastDir = "/usr/local/rmblast/bin"
hmmerDir = "/usr/local/hmmer/bin"
mafftDir = "/usr/local/mafft/bin"
//dialignDir = "/usr/local/dialign-2.2.1"
dialignDir = "/usr/local/dialign-tx-1.0.2"
kalignDir = "/usr/local/kalign2"
tcoffeeDir = "/usr/local/t_coffee/bin"
probconsDir = "/usr/local/probcons/bin"

clustalOmegaDir = "/u1/local/clustal-omega-1.2.4-binary"
fsaDir = "/usr/local/fsa/bin"
muscleDir = "/usr/local/bin"
exonerateDir = "/usr/local/exonerate-2.2.0-x86_64/bin"
repeatmodelerDir = "/home/rhubley/projects/RepeatModeler"

// Generate a series of samples from the orffrags file increasing in size:
//   100, 150, 200, 250, 300, 350, 400, 450, 500, 1000 
seqFile = file( params.familyDir + "/orffrags.fa")
//
curated_protein = file(params.curated_protein)
db_proteins = file(params.db_proteins)

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
log.info "Curated Protein File: " + params.curated_protein
log.info "DB Proteins         : " + params.db_proteins
if ( runMuscle ) {
  log.info "Muscle DIR          : " + muscleDir
}
if ( runMafft ) {
  log.info "Mafft DIR           : " + mafftDir
}
if ( runClustalO ) {
  log.info "ClustalOmega DIR       : " + clustalOmegaDir
}
if ( runTCoffee ) {
  log.info "TCoffee DIR         : " + tcoffeeDir
}
if ( runProbcons ) {
  log.info "ProbCons DIR        : " + probconsDir
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
  file "*sample.fa" into benchmarkFilesForClustalO
  file "*sample.fa" into benchmarkFilesForMafft 
  file "*sample.fa" into benchmarkFilesForMuscle 
  file "*sample.fa" into benchmarkFilesForRefiner 
  file "*sample.fa" into benchmarkFilesForKalign 
  file "*sample.fa" into benchmarkFilesForDialign 
  file "*sample.fa" into benchmarkFilesForFSA 
  file "*sample.fa" into benchmarkFilesForTCoffee
  file "*sample.fa" into benchmarkFilesForProbcons
  file "" into benchmarkFilesForBaseline
  file "*sample.fa"

  script:
  """
  # 
  # Generate samples of various sizes.  
  #   The random seed 1608835988 produced the paper results
  #
  ${workflow.projectDir}/util/sampleFromFA.pl "100,150,200,250,300,350,400,450,500,1000" ${seqFile} 1608835988
  """
}


process runFSA {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForFSA.flatten()

  when:
  runFSA

  output:
  file "*-fsa.fa"
  file "*-fsa.cons.fa"
  file "*-fsa.blastx"
  file "*-fsa.dbprot.blastx"
  file "*-fsa.exonerate"
  file "*-fsa.time"

  script:
  methodPrefix = "fsa"
  sampleSize = (referenceSeqFile.name =~ /^.*\-(\d+)sample.fa/)[0][1]
  """
  if [ ${sampleSize} -lt 350 ] 
  then
    #### Run FSA
    SECONDS=0
    ${fsaDir}/fsa ${referenceSeqFile} > ${referenceSeqFile.baseName}-${methodPrefix}.fa
    echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time
  
    # Gen cons and eval
    ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
    ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  else 
    touch ${referenceSeqFile.baseName}-${methodPrefix}.fa
    touch ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
    echo "0" > ${referenceSeqFile.baseName}-${methodPrefix}.time
    touch ${referenceSeqFile.baseName}-${methodPrefix}.blastx
    touch ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx
    touch ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  fi
  """
}

//  KAlign: https://msa.sbc.su.se/cgi-bin/msa.cgi
//    parameters suggested by website for DNA
process runKalign {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForKalign.flatten()

  when:
  runKalign

  output:
  file "*-kalign.fa"
  file "*-kalign.cons.fa"
  file "*-kalign.blastx"
  file "*-kalign.dbprot.blastx"
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
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}

process runTCoffee {
  cpus 16
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForTCoffee.flatten()

  when:
  runTCoffee

  output:
  file "*-tcoffee.fa"
  file "*-tcoffee.cons.fa"
  file "*-tcoffee.blastx"
  file "*-tcoffee.dbprot.blastx"
  file "*-tcoffee.exonerate"
  file "*-tcoffee.time"

  script:
  methodPrefix = "tcoffee"
  """
  SECONDS=0
  ${tcoffeeDir}/t_coffee -n_core 16 -thread 16 -max_n_proc 16 -seq ${referenceSeqFile} -output fasta_aln -outfile ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}

process runProbcons {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForProbcons.flatten()

  when:
  runProbcons

  output:
  file "*-probcons.fa"
  file "*-probcons.cons.fa"
  file "*-probcons.blastx"
  file "*-probcons.dbprot.blastx"
  file "*-probcons.exonerate"
  file "*-probcons.time"

  script:
  methodPrefix = "probcons"
  """
  SECONDS=0
  ${probconsDir}/probcons ${referenceSeqFile} > ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}


//  DIALIGN: http://dialign.gobics.de/
// Switched to dialigntx because it does appear to perform better
process runDialign {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForDialign.flatten()

  when:
  runDialign

  output:
  file "*-dialign.fa"
  file "*-dialign.cons.fa"
  file "*-dialign.blastx"
  file "*-dialign.dbprot.blastx"
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
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}


process runClustalO {
  cpus 16
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForClustalO.flatten()

  when:
  runClustalO

  output:
  file "*-clustalo.fa"
  file "*-clustalo.cons.fa"
  file "*-clustalo.blastx"
  file "*-clustalo.dbprot.blastx"
  file "*-clustalo.exonerate"
  file "*-clustalo.time"

  script:
  methodPrefix = "clustalo"
  """
  SECONDS=0
  #### Run ClustalOmega
  ${clustalOmegaDir}/clustalo --infile=${referenceSeqFile} --outfile=${referenceSeqFile.baseName}-${methodPrefix}.fa --outfmt=FASTA --threads 16
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}

 

process runMafft {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForMafft.flatten()

  when:
  runMafft

  output:
  file "*-mafft.fa"
  file "*-mafft.cons.fa"
  file "*-mafft.blastx"
  file "*-mafft.dbprot.blastx"
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
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}


process runMuscle {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForMuscle.flatten()

  when:
  runMuscle

  output:
  file "*-muscle.fa"
  file "*-muscle.cons.fa"
  file "*-muscle.blastx"
  file "*-muscle.dbprot.blastx"
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
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
  """
}


process runRefiner {
  publishDir "${outputDir}", mode: 'copy'

  input:
  file curated_protein
  file db_proteins
  file referenceSeqFile from benchmarkFilesForRefiner.flatten()

  when:
  runRefiner

  output:
  file "*-refiner.fa"
  file "*-refiner.cons.fa"
  file "*-refiner.blastx"
  file "*-refiner.dbprot.blastx"
  file "*-refiner.exonerate"
  file "*-refiner.time"

  script:
  methodPrefix = "refiner"
  """
  SECONDS=0
  # Run refiner and generate a filename like gput100-train-refiner.fa for output
  ${repeatmodelerDir}/Refiner ${referenceSeqFile} >& ${referenceSeqFile.baseName}-${methodPrefix}.log
  ## eval Auto Run Blocker
  #${repeatmodelerDir}/util/Linup ${referenceSeqFile}.refiner.stk > alistart
  #${repeatmodelerDir}/util/protocol/AutoRunBlocker.pl -l alistart -w 7 -mc 4 -mr 2 > cons
  #${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  #${repeatmodelerDir}/util/protocol/AutoRunBlocker.pl -l alistart -w 15 -mc 4 -mr 2 > cons
  #${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  #${repeatmodelerDir}/util/protocol/AutoRunBlocker.pl -l alistart -w 24 -mc 4 -mr 2 > cons
  #${repeatmodelerDir}/util/alignAndCallConsensus.pl -c cons -e ${referenceSeqFile}
  #${repeatmodelerDir}/util/protocol/AutoRunBlocker.pl -l alistart -w 5 -mc 4 -mr 2 > cons
  #${repeatmodelerDir}/util/alignAndCallConsensus.pl -re -c cons -e ${referenceSeqFile}
  #${repeatmodelerDir}/util/Linup -msa rep.out > ${referenceSeqFile.baseName}-${methodPrefix}.fa
  ##
  mv ${referenceSeqFile}.refiner_cons ${referenceSeqFile}-refiner.cons.fa
  ${repeatmodelerDir}/util/Linup -msa ${referenceSeqFile}.refiner.stk > ${referenceSeqFile.baseName}-${methodPrefix}.fa
  echo \$SECONDS > ${referenceSeqFile.baseName}-${methodPrefix}.time

  # Gen cons and eval
  ${repeatmodelerDir}/util/Linup -consensus -name ${referenceSeqFile.baseName}-${methodPrefix} ${referenceSeqFile.baseName}-${methodPrefix}.fa > ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa
  ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.blastx 
    ${blastDir}/blastx -query ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa -subject ${db_proteins} > ${referenceSeqFile.baseName}-${methodPrefix}.dbprot.blastx 
  ${exonerateDir}/exonerate ${referenceSeqFile.baseName}-${methodPrefix}.cons.fa ${curated_protein} > ${referenceSeqFile.baseName}-${methodPrefix}.exonerate
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
