##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     sfs.sh                       					###
###   					Estimate sfs using angsd						### 
###########################################################################


#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N sfs


cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load samtools
module load R/3.4.2

# 1D SFS
# only using subsamples #
 
/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
-GL 1 -doPost 3 -nInd 7 -fold 1 -minQ 20 -minMapQ 50 -doMajorMinor 1 -doMaf 1 \
-anc /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1.list \
-pest /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_mqu.sfs \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/az_sub1_mqu

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
-GL 1 -doPost 3 -nInd 7 -fold 1 -minQ 20 -minMapQ 50 -doMajorMinor 1 -doMaf 1 \
-anc /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/tx_s1.list \
-pest /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_mqu.sfs \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/tx_sub1_mqu

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
-GL 1 -doPost 3 -nInd 7 -fold 1 -minQ 20 -minMapQ 50 -doMajorMinor 1 -doMaf 1 \
-anc /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/nm_s1.list \
-pest /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_mqu.sfs \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/nm_sub1_mqu

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/az_sub1_mqu.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/az_sub1_mqu.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/tx_sub1_mqu.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/tx_sub1_mqu.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/nm_sub1_mqu.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/nm_sub1_mqu.sfs



#Bootstrap

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/az_sub1_mqu.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/az_sub1_mqu_boostrap.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/tx_sub1_mqu.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/tx_sub1_mqu_boostrap.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/nm_sub1_mqu.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/pest2/nm_sub1_mqu_boostrap.sfs
