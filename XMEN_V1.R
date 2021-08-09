#############################################
######    Variant Calling Pipeline     ######
#######    Implemented by Min Heo    ########
#######   2021 07-29 ~ 2021-08-08    ########
#############################################

#### Set Env.
setRepositories(ind=1:8)

library(foreach)
library(doParallel)
library(tidyverse)
library(dplyr)
library(parallel)

#### Load file
save.image(file="XMEN.Rdata")
load(file="XMEN.Rdata")

#### Functions 
###############################################################################
##################### 1.Trimming_Trimmomatic Function #########################
###############################################################################

XMEN_Trimmomatic <- function(RealThreads,PutYourTophredNumber,PutYoutFastqDir,id,Full_id_left,Full_id_right,id_left_output_paired,
                             id_right_output_paired,id_left_output_unpaired,id_right_output_unpaired,PutYourTrimmoOutputDir_Pair,
                             PutYourTrimmoOutputDir_UnPair,PutYourTrimAdapter_Dir_ID,PutSeedMismatches,PutPalindromClipThreshold,
                             PutSimpleClipThreshold,PutLEADING,PutTRAILING,PutWindowsize,PutReQuiredQuality,PutMINLEN,Trimming_directory){
  setwd(Trimming_directory)
  TrimmomaticCMD <- paste0("java -jar /data/bin/trimmomatic-0.33.jar PE -threads ",RealThreads," -phred",PutYourTophredNumber," ",PutYoutFastqDir,"/",Full_id_left," ",PutYoutFastqDir,"/",Full_id_right," ",PutYourTrimmoOutputDir_Pair,"/",id_left_output_paired," ",PutYourTrimmoOutputDir_UnPair,"/",id_left_output_unpaired," ",PutYourTrimmoOutputDir_Pair,"/",id_right_output_paired," ",PutYourTrimmoOutputDir_UnPair,"/",id_right_output_unpaired," ","ILLUMINACLIP:",PutYourTrimAdapter_Dir_ID,":",PutSeedMismatches,":",PutPalindromClipThreshold,":",PutSimpleClipThreshold," ","LEADING:",PutLEADING," ","TRAILING:",PutTRAILING," ","SLIDINGWINDOW:",PutWindowsize,":",PutReQuiredQuality," ","MINLEN:",PutMINLEN," ","2> ",Trimming_directory,"/",id,".log")
  system(TrimmomaticCMD)
  print(paste0(id,"'s Trimming process starts"))
}

###############################################################################
########################## 2.QC_fastqC Function ###############################
###############################################################################
XMEN_fastQC <- function(QC_directory,RealThreads,PutYourTrimmoOutputDir_Pair,id_left_output_paired,id_right_output_paired,id){
  setwd(QC_directory)
  FastqcCMD <- paste0("/data/bin/FastQC/fastqc -o ",QC_directory," --extract -f fastq -t ",RealThreads," ",PutYourTrimmoOutputDir_Pair,"/",id_left_output_paired," ",PutYourTrimmoOutputDir_Pair,"/",id_right_output_paired)
  system(FastqcCMD)
}

###############################################################################
########################## 3.Alignment_BWA Function ###########################
###############################################################################
XMEN_BWA <- function(Mapping_directory,RealThreads,PutYourReferenceFile,PutYourTrimmoOutputDir_Pair,id_left_output_paired,id_right_output_paired,id){
  setwd(Mapping_directory)
  bwaCMD <- paste0("/data/bin/bwa mem -t ",RealThreads," -R '@RG\\tID:Korea1k\\tPL:ILLUMINA\\tSM:",id,"\\tLB:WonLab' ",PutYourReferenceFile," ",PutYourTrimmoOutputDir_Pair,"/",id_left_output_paired," ",PutYourTrimmoOutputDir_Pair,"/",id_right_output_paired," > ",Mapping_directory,"/",id,".sam")
  system(bwaCMD)
}

#################################################################################
######################## 4.Converting & Sorting_ Samtools #######################
#################################################################################
XMEN_Samtools_Conv_Sort <- function(Sorting_directory,PutYourReferenceIndexingFile,RealThreads,Mapping_directory,id){
  setwd(Sorting_directory)
  Conv_Sort_CMD <- paste0("/data/bin/samtools view -bt ",PutYourReferenceIndexingFile," -Sb -@ ",RealThreads," ",Mapping_directory,"/",id,".sam | /data/bin/samtools sort -@ ",RealThreads," -o ",Sorting_directory,"/",id,".Sorted.bam")
  system(Conv_Sort_CMD)
}

#################################################################################
############################ 5.Indexing_ Samtools ###############################
#################################################################################
XMEN_Samtools_Indexing <- function(Indexing_directory,Sorting_directory,id,RealThreads){
  setwd(Indexing_directory)
  Indexing_CMD <- paste0("/data/bin/samtools index ",Sorting_directory,"/",id,".Sorted.bam -@ ",RealThreads," -b ",Indexing_directory,"/",id,".index")
  system(Indexing_CMD)
}

#################################################################################
##################### 6.MarkDuplicates_Sortsam_ PICARD ##########################
#################################################################################
XMEN_MarkDup <- function(Markdup_SortSam_directory,Sorting_directory,id){
  setwd(Markdup_SortSam_directory)
  Indexing_CMD <- paste0("java -jar /data/bin/picard/MarkDuplicates.jar I=",Sorting_directory,"/",id,".Sorted.bam O=",Markdup_SortSam_directory,"/",id,".dedup.bam M=",Markdup_SortSam_directory,"/",id,".markdup.metrics.txt")
  system(Indexing_CMD)
}

XMEN_SortSam <- function(Markdup_SortSam_directory,id,RealThreads){
  setwd(Markdup_SortSam_directory)
  SortSam_CMD <- paste0("java -jar /data/bin/picard/SortSam.jar INPUT=",Markdup_SortSam_directory,"/",id,".dedup.bam OUTPUT=",Markdup_SortSam_directory,"/",id,".sortsam.bam SORT_ORDER=coordinate && /data/bin/samtools index ",Markdup_SortSam_directory,"/",id,".sortsam.bam -@ ",RealThreads)
  system(SortSam_CMD)
}

#################################################################################
########################## 7.BQSR_ApplyBQSR_ GATK ###############################
#################################################################################
XMEN_BQSR <- function(BQSR_directory,Markdup_SortSam_directory,id,PutYourReferenceFile,PutYourKnownSitesVCF){
  setwd(BQSR_directory)
  Memory<- system("free -m | grep ^Mem | awk '{print $4}'",intern=T)
  Memory<- as.numeric(Memory)/2
  RealMemory<- Memory%/%fastq_file_length
  BQSR_CMD <- paste0("/data/bin/gatk --java-options '-Xmx",RealMemory,"m -Xms",RealMemory,"m' BaseRecalibrator -I ",Markdup_SortSam_directory,"/",id,".sortsam.bam -R ",PutYourReferenceFile," --known-sites ",PutYourKnownSitesVCF," -O ",BQSR_directory,"/",id,".recal_data.table && /data/bin/gatk --java-options '-Xmx",RealMemory,"m -Xms",RealMemory,"m' ApplyBQSR -I ",Markdup_SortSam_directory,"/",id,".sortsam.bam -R ",PutYourReferenceFile," --bqsr-recal-file ",BQSR_directory,"/",id,".recal_data.table -O ",BQSR_directory,"/",id,".final.bam")
  system(BQSR_CMD)
}
#################################################################################
########################## 8.HaplotypeCaller_ GATK ##############################
#################################################################################
XMEN_HaplotypeCaller <- function(VaraintCalling_directory,PutYourReferenceFile,Markdup_SortSam_directory,id,RealThreads){
  setwd(VaraintCalling_directory)
  Memory<- system("free -m | grep ^Mem | awk '{print $4}'",intern=T)
  Memory<- as.numeric(Memory)/2
  RealMemory<- Memory%/%fastq_file_length
  HaplotypeCaller_CMD <- paste0("/data/bin/gatk --java-options '-Xmx",RealMemory,"m -Xms",RealMemory,"m' HaplotypeCaller -R ",PutYourReferenceFile," -I ",BQSR_directory,"/",id,".final.bam -O ",VaraintCalling_directory,"/",id,".vcf -ERC GVCF --output-mode EMIT_VARIANTS_ONLY -ploidy 2 --stand-call-conf 2 --native-pair-hmm-threads ",RealThreads)
  system(HaplotypeCaller_CMD)
}

#################################################################################
######################## 9.Making Thread_function ###############################
#################################################################################
realthread <- function(PutYourThreadsNumber,fastq_file_length){
  realnumber <- PutYourThreadsNumber %/% fastq_file_length
  return(realnumber)
}

#### Input Indexes
Args <- commandArgs(trailingOnly = T)
#### Necessary idx (not have default value)
idxVCF_pipeline_Directory <- which(Args=="-Dp")
idxPutYoutFastqDir <- which(Args=="-Df")
idxPutYourThreadsNumber <- which(Args=="-t")
idxPutYourTrimAdapter_Dir_ID <- which(Args=="-a")
idxPutYourReferenceFile <- which(Args=="-R")
idxPutYourReferenceIndexingFile <- which(Args=="-Ri")
idxPutYourKnownSitesVCF <- which(Args=="-K")
#### Unnecessary idx (have default value)
idxTophred <- which(Args=="--p")
idxMismatches <- which(Args=="--m")
idxPalinClipThresh <- which(Args=="--pthr")
idxSimpleClipThresh <- which(Args=="--sthr")
idxLEADING <- which(Args=="--l")
idxTRAILING <- which(Args=="--t")
idxWindowsize <- which(Args=="--w")
idxReQuiredQuality <- which(Args=="--q")
idxMINLEN <- which(Args=="--ml")

#### Input file
#### Necessary (not have default value)
VCF_pipeline_Directory <- Args[idxVCF_pipeline_Directory+1] # "/home2/mheo/VariantCalling_Pipeline/PipelineTest_HM"  
PutYoutFastqDir <- Args[idxPutYoutFastqDir+1] # "/home2/mheo/VariantCalling_Pipeline/PipelineTest_HM/Toydata"
PutYourThreadsNumber <- as.numeric(Args[idxPutYourThreadsNumber+1]) #12 
PutYourTrimAdapter_Dir_ID <- Args[idxPutYourTrimAdapter_Dir_ID+1] # "/data/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa" 
PutYourReferenceFile <- Args[idxPutYourReferenceFile+1] # "/home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/ucsc.hg19.chr.only.fasta" 
PutYourReferenceIndexingFile <- Args[idxPutYourReferenceIndexingFile+1] # "/home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/ucsc.hg19.chr.only.fasta.fai" 
PutYourKnownSitesVCF <- Args[idxPutYourKnownSitesVCF+1] # "/home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/dbsnp_151.hg19.vcf" 
#### default value
PutYourTophredNumber <- 33
PutSeedMismatches <- 2
PutPalindromClipThreshold <- 30
PutSimpleClipThreshold <- 10
PutLEADING <- 3
PutTRAILING <- 3
PutWindowsize <- 4
PutReQuiredQuality <- 15
PutMINLEN <- 36

WORK_DIRECTORY <- VCF_pipeline_Directory
setwd(WORK_DIRECTORY)

#### If you input specific value(arguments) for your study design, these values are change. (Unnecessary Values)
if(length(idxTophred)>0){
  PutYourTophredNumber <- Args[idxTophred+1]
}
if(length(idxMismatches)>0){
  PutSeedMismatches <- Args[idxMismatches+1] 
}
if(length(idxPalinClipThresh)>0){
  PutPalindromClipThreshold <- Args[idxPalinClipThresh+1] 
}
if(length(idxSimpleClipThresh)>0){
  PutSimpleClipThreshold <- Args[idxSimpleClipThresh+1]
}
if(length(idxLEADING)>0){
  PutLEADING <- Args[idxLEADING+1] 
}
if(length(idxTRAILING)>0){
  PutTRAILING <- Args[idxTRAILING+1]
}
if(length(idxWindowsize)>0){
  PutWindowsize <- Args[idxWindowsize+1] 
}
if(length(idxReQuiredQuality)>0){
  PutReQuiredQuality <- Args[idxReQuiredQuality+1] 
}
if(length(idxMINLEN)>0){
  PutMINLEN <- Args[idxMINLEN+1]
}

### Make Directory 
MKdir <- c()
for(i in 1:10){
  string <- if(i == 1){"/1.Trimmming_Trimmomatic"}
  else if(i == 2){"/2.QC_FastQC"}
  else if(i == 3){"/3.Alignment_BWA"}
  else if(i == 4){"/4.Converting_Sorting_Samtools"}
  else if(i == 5){"/5.BAM_Indexing_Samtools"}
  else if(i == 6){"/6.MarkDup_SortSam_GATK"}
  else if(i == 7){"/7.BQSR_GATK"}
  else if(i == 8){"/8.VariantCalling_GATK"}
  else if(i == 9){"/1.Trimmming_Trimmomatic/Paired"}
  else if(i == 10){"/1.Trimmming_Trimmomatic/Unpaired"}
  MKdir[i] <- paste0("mkdir ",VCF_pipeline_Directory,string)
}

for(i in 1:length(MKdir)){
  system(MKdir[i])
}
### Program file 
Trimming_directory <- MKdir[1] %>% str_remove("mkdir ")
PutYourTrimmoOutputDir_Pair <- MKdir[9] %>% str_remove("mkdir ")
PutYourTrimmoOutputDir_UnPair <- MKdir[10] %>% str_remove("mkdir ")
QC_directory <- MKdir[2] %>% str_remove("mkdir ")
Mapping_directory <- MKdir[3] %>% str_remove("mkdir ")
Sorting_directory <- MKdir[4] %>% str_remove("mkdir ")
Indexing_directory <- MKdir[5] %>% str_remove("mkdir ")
Markdup_SortSam_directory <- MKdir[6] %>% str_remove("mkdir ")
BQSR_directory <- MKdir[7] %>% str_remove("mkdir ")
VaraintCalling_directory <- MKdir[8] %>% str_remove("mkdir ")

#######################################################################
#####################    ID Making Process   ##########################
#######################################################################
FastqFile <- list.files(path=PutYoutFastqDir,pattern="R1.clean.fq.gz")
fastq_file_length <- length(FastqFile)

id_r <- strsplit(FastqFile,"_")
id_unlist <- as.data.frame(id_r)

id <- c()
id_left <- c()
id_right<- c()
Full_id_left <- c()
Full_id_right <- c()
id_left_output_paired<- c()
id_right_output_paired<- c()
id_left_output_unpaired<- c()
id_right_output_unpaired<- c()

for(i in 1:fastq_file_length){
  id[i] <- as.character(paste0(id_unlist[1,i],"_",id_unlist[2,i]))
  
  id_left[i] <- paste0(id[i],"_R1")
  id_right[i] <- paste0(id[i],"_R2")
  
  Full_id_left[i] <- paste0(id[i],"_R1.clean.fq.gz")
  Full_id_right[i] <- paste0(id[i],"_R2.clean.fq.gz")
  
  id_left_output_paired[i] <-  paste0(id_left[i],".paired.fq.gz")  
  id_right_output_paired[i] <- paste0(id_right[i],".paired.fq.gz") 
  
  id_left_output_unpaired[i] <-  paste0(id_left[i],".unpaired.fq.gz")  
  id_right_output_unpaired[i] <-  paste0(id_right[i],".unpaired.fq.gz")
}

#################################################################################
############################# Execution Part ####################################
#################################################################################

#### Parallel Setting
Cl <- makeCluster(PutYourThreadsNumber/2) 
doParallel::registerDoParallel(Cl)

#### Making Real Thread Numbers 
RealThreads <- realthread(PutYourThreadsNumber,fastq_file_length)

#### Execution part 
foreach(i=1:fastq_file_length,.packages = c("tidyverse","dplyr")) %dopar% {
  ## Trimmomatic (Trimming) 
  if(length(list.files(path=PutYourTrimmoOutputDir_Pair,pattern=paste0(id_left[i],".paired.fq.gz")))!=0){
    system("echo 'Data is existed ... Skip Trimming process'")
  }else{
    if(length(list.files(path=PutYoutFastqDir,pattern=paste0(id_left[i],".clean.fq.gz")))!=0 && length(PutYourTrimAdapter_Dir_ID)!= 0){
      system("echo 'Start Trimmomatic'")
      XMEN_Trimmomatic(RealThreads,PutYourTophredNumber,PutYoutFastqDir,id[i],Full_id_left[i],Full_id_right[i],id_left_output_paired[i],
                       id_right_output_paired[i],id_left_output_unpaired[i],id_right_output_unpaired[i],PutYourTrimmoOutputDir_Pair,PutYourTrimmoOutputDir_UnPair,
                       PutYourTrimAdapter_Dir_ID,PutSeedMismatches,PutPalindromClipThreshold,PutSimpleClipThreshold,PutLEADING,PutTRAILING,
                       PutWindowsize,PutReQuiredQuality,PutMINLEN,Trimming_directory)
    }else{
      system("echo 'Error in your inputs for processing Trimmomatic'")
      system("echo 'Please check your inputs & directories & files related to trimming'")
      stop()
    }
  }
  ## Fastqc (QC checking)
  if(length(list.files(path=QC_directory,pattern=paste0(id_left[i],".paired_fastqc.html")))!=0){
    system("echo 'Data is existed ... Skip fastQC process'")
  }else{
    if(length(list.files(path=PutYourTrimmoOutputDir_Pair,pattern=paste0(id_left[i],".paired.fq.gz")))!=0){
      system("echo 'Start fastQC'")
      XMEN_fastQC(QC_directory,RealThreads,PutYourTrimmoOutputDir_Pair,id_left_output_paired[i],id_right_output_paired[i],id[i])
    }else{
      system("echo 'Error in your inputs for processing fastQC'")
      system("echo 'Please check your inputs & directories & files related to fastQC'")
      stop()
    }
  }
  ## BWA (Mapping)
  if(length(list.files(path=Mapping_directory,pattern=paste0(id[i],".sam")))!=0){
    system("echo 'Data is existed ... Skip BWA process'")
  }else{
    if(length(list.files(path=PutYourTrimmoOutputDir_Pair,pattern=paste0(id_left[i],".paired.fq.gz")))!=0 && length(PutYourReferenceFile)!=0){
      system("echo 'Start BWA'")
      XMEN_BWA(Mapping_directory,RealThreads,PutYourReferenceFile,PutYourTrimmoOutputDir_Pair,id_left_output_paired[i],id_right_output_paired[i],id[i])
    }else{
      system("echo 'Error in your inputs for processing BWA'")
      system("echo 'Please check your inputs & directories & files related to BWA'")
      stop()
    }
  }
  ## Samtools (Converting & Sorting)
  if(length(list.files(path=Sorting_directory,pattern=paste0(id[i],".Sorted.bam")))!=0){
    system("echo 'Data is existed ... Skip Convertinng & Sorting process'")
  }else{
    if(length(list.files(path=Mapping_directory,pattern=paste0(id[i],".sam")))!=0){
      system("echo 'Start Converting & Sorting'")
      XMEN_Samtools_Conv_Sort(Sorting_directory,PutYourReferenceIndexingFile,RealThreads,Mapping_directory,id[i])
    }else{
      system("echo 'Error in your inputs for processing Converting & Sorting'")
      system("echo 'Please check your inputs & directories & files related to Converting & Sorting'")
      stop()
    }
  }
  ## Samtools Indexing (BAM file Indexing)
  if(length(list.files(path=Indexing_directory,pattern=paste0(id[i],".index")))!=0){
    system("echo 'Data is existed ... Skip BAM file Indexing process'")
  }else{
    if(length(list.files(path=Sorting_directory,pattern=paste0(id[i],".Sorted.bam")))!=0){
      system("echo 'Start BAM Indexing'")
      XMEN_Samtools_Indexing(Indexing_directory,Sorting_directory,id[i],RealThreads)
    }else{
      system("echo 'Error in your inputs for processing BAM file Indexing'")
      system("echo 'Please check your inputs & directories & files related to BAM file Indexing'")
      stop()
    }
  }
  ## Picard MarkDuplicates (Remove Duplicates)
  if(length(list.files(path=Markdup_SortSam_directory,pattern=paste0(id[i],".dedup.bam")))!=0 && length(list.files(path=Markdup_SortSam_directory,pattern=".markdup.metrics.txt"))!=0){
    system("echo 'Data is existed ... Skip MarkDuplicates process'")
  }else{
    if(length(list.files(path=Sorting_directory,pattern=paste0(id[i],".Sorted.bam")))!=0){
      system("echo 'Start MarkDuplicates'")
      XMEN_MarkDup(Markdup_SortSam_directory,Sorting_directory,id[i])
    }else{
      system("echo 'Error in your inputs for processing MarkDuplicates'")
      system("echo 'Please check your inputs & directories & files related to MarkDuplicates'")
      stop()
    }
  }
  ## Picard SortSam (Sorting dedup.bam file)
  if(length(list.files(path=Markdup_SortSam_directory,pattern=paste0(id[i],".sortsam.bam")))!=0){
    system("echo 'Data is existed ... Skip SortSam process'")
  }else{
    if(length(list.files(path=Markdup_SortSam_directory,pattern=paste0(id[i],".dedup.bam")))!=0 && length(list.files(path=Markdup_SortSam_directory,pattern=".markdup.metrics.txt"))!=0){
      system("echo 'Start Sortsam'")
      XMEN_SortSam(Markdup_SortSam_directory,id[i],RealThreads)
    }else{
      system("echo 'Error in your inputs for processing SortSam'")
      system("echo 'Please check your inputs & directories & files related to SortSam'")
      stop()
    }
  }
  ## GATK BaseRecalibrator(BQSR) & Apply BQSR (Base Recalibration process)
  if(length(list.files(path=BQSR_directory,pattern=paste0(id[i],".final.bam")))!=0){
    system("echo 'Data is existed ... Skip BQSR process'")
  }else{
    if(length(list.files(path=Markdup_SortSam_directory,pattern=paste0(id[i],".sortsam.bam")))!=0 && length(PutYourReferenceFile)!=0 && length(PutYourKnownSitesVCF)!=0){
      system("echo 'Start BQSR'")
      XMEN_BQSR(BQSR_directory,Markdup_SortSam_directory,id[i],PutYourReferenceFile,PutYourKnownSitesVCF)
    }else{
      system("echo 'Error in your inputs for processing BQSR'")
      system("echo 'Please check your inputs & directories & files related to BQSR'")
      stop()
    }
  }
  ## GATK HaplotypeCaller (BAM -> gVCF)
  if(length(list.files(path=VaraintCalling_directory,pattern=paste0(id[i],".vcf")))!=0){
    system("echo 'Data is existed ... Skip HaplotypeCaller process'")
  }else{
    if(length(list.files(path=BQSR_directory,pattern=paste0(id[i],".final.bam")))!=0 && length(PutYourReferenceFile)!=0){
      system("echo 'Start HaplotypeCaller'")
      XMEN_HaplotypeCaller(VaraintCalling_directory,PutYourReferenceFile,Markdup_SortSam_directory,id[i],RealThreads)
      if(i==fastq_file_length){
        system("echo 'All processes are over'")
      }
    }else{
      system("echo 'Error in your inputs for processing HaplotypeCaller'")
      system("echo 'Please check your inputs & directories & files related to HaplotypeCaller'")
      stop()
    }
  }
}
parallel:::stopCluster(Cl)


