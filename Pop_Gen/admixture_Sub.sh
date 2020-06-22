##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     admixture_Sub.sh                       			###
###   				admixture analysis using GL 	 					### 
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=100:00:00
#PBS -N admix_sub

cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load samtools

j="1"
while [ $j -lt 11 ]
do
	i="2"
	while [ $i -lt 11 ]
	do
		/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/NGSadmix -P 80 \
		-K ${i} -minMaf 0.05 -maxiter 50000 -tol 1e-9 -tolLike50 1e-9 \
		-likes /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/sub1.beagle.gz \
		-outfiles /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/admix/r${j}/sub1_k${i} 
 		i=$[$i+1]
	done
	j=$[$j+1]
done