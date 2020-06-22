##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     genotyple_likelihood_sub.sh                     ###
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N sub_gl

cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load samtools
module load R/3.4.2

# GL #

#all_pop
/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
-nInd 21 -minQ 20 -minMapQ 50 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 -doDepth 1 -doCounts 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 15 \
-SNP_pval 1e-6 -setMaxDepth 100 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/subsample1.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/depth/sub1

#by_pop

cd /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd sites index \
/scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/sub1.sites


/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-nInd 7 -minQ 20 -minMapQ 50 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 5 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-sites /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/sub1.sites \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/az_sub1

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-nInd 7 -minQ 20 -minMapQ 50 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 5 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-sites /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/sub1.sites \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/tx_s1.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/tx_sub1

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-nInd 7 -minQ 20 -minMapQ 50 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 5 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-sites /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/sub1.sites \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/nm_s1.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/nm_sub1