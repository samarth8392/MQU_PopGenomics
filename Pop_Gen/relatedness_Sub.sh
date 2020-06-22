##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     relatedness_sub.sh                     			###
###   					Estimate relatedness using ibsrelate 			### 
###########################################################################

#!/bin/sh -l
#PBS -q fnrdewoody
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N relatedness_sub


cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load samtools
module load R/3.4.2

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doGlf 1 -nInd 21 -doHWE 1 \
-doMajorMinor 1 -doMaf 3 \
-minMapQ 50 -minQ 20 -uniqueOnly 1 \
-minHWEpval 0.1 \
-minInd 15 -setMinDepthInd 1 \
-SNP_pval 1e-6 -setMaxDepth 100 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/subsample1.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/relate/sub1

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/ibs \
-glf /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/relate/sub1.glf.gz \
-model 0 \
-nInd 21 -allpairs 1 \
-outFileName /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/relate/sub1

Rscript \
  -e "source('/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/R/read_IBS.R')" \
  -e "res = do_derived_stats(read_ibspair_model0('/scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/relate/sub1.ibspair'))" \
  -e "print(res[,c('ind1', 'ind2', 'nSites', 'Kin', 'R0', 'R1') ])" \
  > /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/relate/sub1_relate.txt