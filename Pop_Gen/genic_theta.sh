##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     genic_theta.sh                       			###
### Estimating nucleotide diversity for MQU genes						###
###########################################################################


#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N genic_theta


cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load samtools
module load R/3.4.2

# 1D SFS

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -nInd 7 -fold 1 -minQ 20 -minMapQ 50 \
-anc /scratch/snyder/m/mathur20/MQU/2020/adapt/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
-bam /scratch/snyder/m/mathur20/MQU/2020/adapt/az_s1.list \
-out /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -nInd 7 -fold 1 -minQ 20 -minMapQ 50 \
-anc /scratch/snyder/m/mathur20/MQU/2020/adapt/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
-bam /scratch/snyder/m/mathur20/MQU/2020/adapt/tx_s1.list \
-out /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -nInd 7 -fold 1 -minQ 20 -minMapQ 50 \
-anc /scratch/snyder/m/mathur20/MQU/2020/adapt/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
-bam /scratch/snyder/m/mathur20/MQU/2020/adapt/nm_s1.list \
-out /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_sub1_genic

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_sub1_genic.saf.idx -maxIter 10000 \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_sub1_genic.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 100 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_genic_mqu_boostrap.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 100 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_genic_mqu_boostrap.sfs

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 100 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_sub1_genic.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_genic_mqu_boostrap.sfs

# Theta #

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -doThetas 1 -nind 7 -fold 1 \
-pest /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic.sfs \
-bam /scratch/snyder/m/mathur20/MQU/2020/adapt/az_s1.list \
-anc /scratch/snyder/m/mathur20/MQU/2020/adapt/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
-out /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/az_sub1_genic

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -doThetas 1 -nind 7 -fold 1 \
-pest /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic.sfs \
-bam /scratch/snyder/m/mathur20/MQU/2020/adapt/tx_s1.list \
-anc /scratch/snyder/m/mathur20/MQU/2020/adapt/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
-out /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/tx_sub1_genic

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 80 \
-GL 1 -doSaf 1 -doThetas 1 -nind 7 -fold 1 \
-pest /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_sub1_genic.sfs \
-bam /scratch/snyder/m/mathur20/MQU/2020/adapt/nm_s1.list \
-anc /scratch/snyder/m/mathur20/MQU/2020/adapt/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
-out /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/nm_sub1_genic

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat print \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/az_sub1_genic.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/az_sub1_genic_persite.txt

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat print \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/tx_sub1_genic.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/tx_sub1_genic_persite.txt

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat print \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/nm_sub1_genic.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/nm_sub1_genic_persite.txt

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/az_sub1_genic.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/az_sub1_genic.txt

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/tx_sub1_genic.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/tx_sub1_genic.txt

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/nm_sub1_genic.thetas.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/theta/nm_sub1_genic.txt


