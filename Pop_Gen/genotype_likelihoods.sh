###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     genotype_likelihoods.sh                       	###
###   					Estimate genotype likelihoods					### 
###########################################################################

#!/bin/sh -l
#PBS -q fnrgenetics
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N gl_all_mqu_depth

cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load samtools


#all_pop
/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
-nInd 74 -minQ 20 -minMapQ 50 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 -doDepth 1 -doCounts 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 60 \
-SNP_pval 1e-6 -setMaxDepth 500 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/bam/bamlist/all_mqu_wgs.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/depth/all_pop_mqu

#get polymorhpic sites from all_pop

cd /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd sites index \
/scratch/snyder/m/mathur20/MQU/2020/angsd/gl/all_pop_chick.sites

#by_pop (Using polymorphic sites from all_pop)

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-nInd 52 -minQ 20 -minMapQ 50 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 1 \
-SNP_pval 1e-6 -setMaxDepth 500 \
-minMaf 0.05 -minHWEpval 0.01 \
-sites /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/all_pop_chick.sites \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/bam/bamlist/az_chick_wgs.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/az_pop_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-nInd 15 -minQ 20 -minMapQ 50 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 1 \
-SNP_pval 1e-6 -setMaxDepth 500 \
-minMaf 0.05 -minHWEpval 0.01 \
-sites /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/all_pop_chick.sites \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/bam/bamlist/tx_chick_wgs.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/tx_pop_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-nInd 7 -minQ 20 -minMapQ 50 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 1 \
-SNP_pval 1e-6 -setMaxDepth 500 \
-minMaf 0.05 -minHWEpval 0.01 \
-sites /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/all_pop_chick.sites \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/bam/bamlist/nm_chick_wgs.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/nm_pop_chick