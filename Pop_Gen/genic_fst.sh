##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     genic_fst.sh                       				###
### Estimating global FST between pops in MQU genes						###
###########################################################################


#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N genic_fst


cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load samtools
module load R/3.4.2

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 -bootstrap 100 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic.saf.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/az.tx.genic.ml

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 -bootstrap 100 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_sub1_genic.saf.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/tx.nm.genic.ml

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 -bootstrap 100 \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic.saf.idx \
/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/nm_sub1_genic.saf.idx \
> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/az.nm.genic.ml

# Calculate global Fst for each bootrap 2D-SFS

#AZ-TX (Do it similarly for TX-NM and AZ-NM) #

i="1"
while read line
do
	echo ${line[0]} > /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/aztx/az_tx_sfs_${i}
	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst index -P 80 \
	/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/az_sub1_genic.saf.idx \
	/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/sfs/tx_sub1_genic.saf.idx \
	-sfs  /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/aztx/az_tx_sfs_${i} \
	-fstout /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/aztx/az_tx_genic_${i}

	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS fst stats \
	/scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/aztx/az_tx_genic_${i}.fst.idx \
	> /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/aztx/az_tx_genic_${i}_globalFst
	
	i=$[$i+1]
done < /scratch/snyder/m/mathur20/MQU/2020/adapt/angsd/fst/bootstrap/az.tx.genic.ml


