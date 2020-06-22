##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     pbs.sh                       					###
###   		Genome wide Population Branch Statistics analysis			### 
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N pbs_genic

cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load samtools
module load R/3.6.1

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst index -P 80 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_sub1_genic.saf.idx \
-sfs /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/az.tx.genic.ml \
-sfs /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/tx.nm.genic.ml \
-sfs /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/az.nm.genic.ml \
-fstout /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/aztxnm_genic

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats \
/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/aztxnm_genic.fst.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/aztxnm_genic_globalpbs


/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats2 \
/scratch/snyder/m/mathur20/MQU/2020/angsd/fst/aztxnm_genic.fst.idx \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/fst/aztxnm_genic.pbs
