# WGS Variant Calling Pipeline by R language

#### Manual for XMEN made by Min Heo.

## Pipeline Flow 

##### Trimmomatic -> fastQC -> BWA -> Samtools -> Picard -> GATK

#### Inputs' index, if index is "-a", it is essential index, and if index is "--a", it is not essential, it has default value. If you want to change these values, you can change them.

## Necessary idx (not have default value, You should enter these values.)
-Dp: Put your Variant Calling Pipeline directory which you want to store results.

-Df: Put your directory that includes fastq files(sample files).

-t: Put your threads number.

-a: Put your adpater file for Trimmomatic.

-R: Put your Reference file. 

-Ri: Put your Reference indexing file.

-K: Puy your known sites for base recalibration.

## Unnecessary idx (have default value,Values in () is default values for implementation , You don't need to enter these values.)

--p: Tophred number (33)

--m: Trimmomatic's mismatch numbers (2)

--pthr: Trimmomatic's Palindrome Clip Threshold (30)

--sthr: Trimmomatic's Simple Clip Threshold (10)

--l: Trimmomatic's Leading (3)

--t: Trimmomatic's Trailing (3)

--w: Trimmomatic's Window Sliding Size (4)

--q: Trimmomatic's requiredQuality (15)

--ml: Trimmomatic's MINLENGTH (36)

## Execution Code Example (Toy Data)
time Rscript /home2/mheo/VariantCalling_Pipeline/PipelineTest_HM/XMEN_V1.5_NoQC.R -Dp /home2/mheo/VariantCalling_Pipeline/hi -Df /home2/mheo/VariantCalling_Pipeline/PipelineTest_HM/Toydata -t 60 -a /data/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa -R /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/ucsc.hg19.chr.only.fasta -Ri /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/ucsc.hg19.chr.only.fasta.fai -K /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/dbsnp_151.hg19.vcf

## Execution Code Example (Real Data)
time Rscript /home2/mheo/VariantCalling_Pipeline/PipelineTest_HM/XMEN_V1.5_NoQC.R -Dp /data6/20210809_VariantCalling/3.Phase3/ -Df /data6/20210809_VariantCalling/0.Rawfile/3.Phase3/ -t 60 -a /data/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa -R /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/ucsc.hg19.chr.only.fasta -Ri /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/ucsc.hg19.chr.only.fasta.fai -K /home2/mheo/VariantCalling_Pipeline/Korea1K_VCF_example/0.Reference/hg19/dbsnp_151.hg19.vcf
