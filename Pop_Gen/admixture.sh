##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     admixture.sh                       				###
###   				admixture analysis using GL 	 					### 
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N admix

cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load samtools

j="1"
while [ $j -lt 11 ]
do
	for i in 1 11
	do
		/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/NGSadmix -P 80 \
		-K ${i} -minMaf 0.05 -maxiter 50000 -tol 1e-9 -tolLike50 1e-9 \
		-likes /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/all_pop_mqu.beagle.gz \
		-outfiles /scratch/snyder/m/mathur20/MQU/2020/angsd/admix/r${j}/all_pop_mqu_k${i} 
 		i=$[$i+1]
	done
	j=$[$j+1]
done
