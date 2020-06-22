
##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                    genotype_calling.sh                       		###
###   			Calling SNP genotypes from whole genome data 			### 
###########################################################################


#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N snp_call

cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load samtools
module load bedops
module load vcftools

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-doBcf 1 -GL 1 -doCounts 1 -nInd 7 -doHWE 1 -doPost 1 -doMajorMinor 1 -doMaf 3 \
--ignore-RG 0 -dogeno 1 -geno_minDepth 5 -minMapQ 50 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-setMinDepthInd 5 -minInd 1 \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1_mqu.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-doBcf 1 -GL 1 -doCounts 1 -nInd 7 -doHWE 1 -doPost 1 -doMajorMinor 1 -doMaf 3 \
--ignore-RG 0 -dogeno 1 -geno_minDepth 5 -minMapQ 50 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-setMinDepthInd 5 -minInd 1 \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1_mqu.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/tx_s1

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-doBcf 1 -GL 1 -doCounts 1 -nInd 7 -doHWE 1 -doPost 1 -doMajorMinor 1 -doMaf 3 \
--ignore-RG 0 -dogeno 1 -geno_minDepth 5 -minMapQ 50 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-setMinDepthInd 5 -minInd 1 \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1_mqu.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/nm_s1

bcftools convert -O v -o /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1.vcf \
/scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1.bcf

bcftools convert -O v -o /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/tx_s1.vcf \
/scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/tx_s1.bcf

bcftools convert -O v -o /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/nm_s1.vcf \
/scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/nm_s1.bcf