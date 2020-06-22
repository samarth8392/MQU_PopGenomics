##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     fst.sh                       					###
###   					Genome wide Fst analysis						### 
###########################################################################


#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N fst


cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load samtools
module load R/3.4.2

# sub Samples #
# 1D SFS using chicken #
/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -nInd 7 -minQ 20 -minMapQ 50 \
-anc /scratch/snyder/m/mathur20/MQU/2020/align/chicken_genome.fa \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1_chick.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -nInd 7 -minQ 20 -minMapQ 50 \
-anc /scratch/snyder/m/mathur20/MQU/2020/align/chicken_genome.fa \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/tx_s1_chick.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -nInd 7 -minQ 20 -minMapQ 50 \
-anc /scratch/snyder/m/mathur20/MQU/2020/align/chicken_genome.fa \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/nm_s1_chick.list \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_chick.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_chick.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_chick.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_chick.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_chick.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_chick.sfs

#calculate the 2dsfs prior

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 -bootstrap 100 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_chick.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_chick.saf.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/az.tx.chick.ml

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 -bootstrap 100 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_chick.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_chick.saf.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/txnm/tx.nm.chick.ml

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 -bootstrap 100 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_chick.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_chick.saf.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/aznm/az.nm.chick.ml

# Calculate global Fst for each bootrap 2D-SFS

#AZ-TX (Do it similarly for TX-NM and AZ-NM)
i="1"
while read line
do
	echo ${line[0]} > /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/txnm/az_tx_sfs_${i}
	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst index -P 200 \
	/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_chick.saf.idx \
	/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_chick.saf.idx \
	-sfs  /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/aztx/az_tx_sfs_${i} \
	-fstout /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/aztx/az_tx_s1_${i}

	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats \
	/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/aztx/az_tx_s1_${i}.fst.idx \
	> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/aztx/az_tx_genic_${i}_globalFst
	
	i=$[$i+1]
done < /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/bootstrap/az.tx.chick.ml

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst index -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_chick.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_chick.saf.idx \
-sfs /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az.tx.chick.ml \
-fstout /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_tx_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst index -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_chick.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_chick.saf.idx \
-sfs /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/tx.nm.chick.ml \
-fstout /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/tx_nm_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst index -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_chick.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_chick.saf.idx \
-sfs /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az.nm.chick.ml \
-fstout /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_nm_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats \
/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_tx_chick.fst.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_tx_chick_globalFst

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats \
/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/tx_nm_chick.fst.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/tx_nm_chick_globalFst

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats \
/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_nm_chick.fst.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_nm_chick_globalFst

#Fst sliding window (window size = 50kb, step = 10kb)
/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats2 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_tx_chick.fst.idx -win 100000 -step 50000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_tx_chick.fst_10050 

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats2 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/tx_nm_chick.fst.idx -win 100000 -step 50000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/tx_nm_chick.fst_10050  

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats2 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_nm_chick.fst.idx -win 100000 -step 50000 \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/az_nm_chick.fst_10050 

