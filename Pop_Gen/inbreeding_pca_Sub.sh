##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     inbreeding_sub.sh                        		###
###   			Individual inbreeding coefficients						### 
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=5,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N inbreeding

cd $PBS_O_WORKDIR

module load use.own
module load conda-env/pcangsd-py2.7.15

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -inbreed 1 -inbreedSites -iter 5000 -maf_iter 5000 -inbreed_iter 5000 -inbreed_tole 1e-9 \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/sub1.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/all_sub1 \
&> /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/all_sub1.log

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -inbreed 1 -inbreedSites -iter 5000 -maf_iter 5000 -inbreed_iter 5000 -inbreed_tole 1e-9 \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/az_sub1.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/az_sub1 \
&> /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/az_sub1.log

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -inbreed 1 -inbreedSites -iter 5000 -maf_iter 5000 -inbreed_iter 5000 -inbreed_tole 1e-9 \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/tx_sub1.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/tx_sub1 \
&> /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/tx_sub1.log

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -inbreed 1 -inbreedSites -iter 5000 -maf_iter 5000 -inbreed_iter 5000 -inbreed_tole 1e-9 \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/nm_sub1.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/nm_sub1 \
&> /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/nm_sub1.log

#PCA

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -iter 5000 -maf_iter 5000 \
-hwe /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/all_sub1.lrt.sites.npy \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/gl/sub1.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/inbreeding/all_sub1_pca