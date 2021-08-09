## WGS Variant Calling Pipeline by R language


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



#### Execution Code Example
Rscript /home2/mheo/VariantCalling_Pipeline/PipelineTest_HM/XMEN_V1.R -Dp /home2/mheipelineTest_HM/XMEN.R -Dp /home2/mheo/VariantCalling_Pipeline/hi/ -Df /home2/mhe o/VariantCalling_Pipeline/PipelineTest_HM/Toydata/ -t 12 -a /data/Trimmomatic-0.933/adapters/TruSeq3-PE-2.fa -R /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/ucsc.hg19.chr.only.fasta -Ri /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/ucsc.hg19.chr.only.fasta.fai -K /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/dbsnp_151.hg19.vcf
