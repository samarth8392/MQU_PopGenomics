##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     theta.sh    									###
###   					Estimate theta from sfs using angsd				### 
###########################################################################


#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N theta


cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load samtools
module load R/3.4.2


/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
-GL 1 -doSaf 1 -doThetas 1 -nind 7 \
-pest /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/az_sub1_mqu.sfs \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1_mqu.list \
-anc /scratch/snyder/m/mathur20/MQU/2020/align/mquen_genome.fa \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/az_s1_mqu

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
-GL 1 -doSaf 1 -doThetas 1 -nind 7 \
-pest /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/tx_sub1_mqu.sfs \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/tx_s1_mqu.list \
-anc /scratch/snyder/m/mathur20/MQU/2020/align/mquen_genome.fa \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/tx_s1_mqu

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
-GL 1 -doSaf 1 -doThetas 1 -nind 7 \
-pest /scratch/snyder/m/mathur20/MQU/2020/angsd/sfs/nm_sub1_mqu.sfs \
-bam /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/nm_s1_mqu.list \
-anc /scratch/snyder/m/mathur20/MQU/2020/align/mquen_genome.fa \
-out /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/nm_s1_mqu

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/az_s1_mqu.thetas.idx \
-win 100000 -step 50000 \
-outnames /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/az_s1_mqu_10050

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/az_s1_mqu.thetas.idx \
-outnames /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/az_s1_mqu

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/tx_s1_mqu.thetas.idx \
-win 100000 -step 50000 \
-outnames /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/tx_s1_mqu_10050

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/tx_s1_mqu.thetas.idx \
-outnames /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/tx_s1_mqu


/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/nm_s1_mqu.thetas.idx \
-win 100000 -step 50000 \
-outnames /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/nm_s1_mqu_10050

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/nm_s1_mqu.thetas.idx \
-outnames /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/nm_s1_mqu


/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat print \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/az_s1_mqu.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/az_s1_mqu_persite.txt

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat print \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/tx_s1_mqu.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/tx_s1_mqu_persite.txt

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat print \
/scratch/snyder/m/mathur20/MQU/2020/angsd/theta/nm_s1_mqu.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/theta/nm_s1_mqu_persite.txt
