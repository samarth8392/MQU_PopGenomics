##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     inbreeding.sh                        			###
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

# inbreeding #

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -inbreed 1 -inbreedSites -iter 5000 -maf_iter 5000 -inbreed_iter 5000 -inbreed_tole 1e-9 \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/all_pop_mqu.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/all_pop_mqu \
&> /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/all_pop_mqu.log

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -inbreed 1 -inbreedSites -iter 5000 -maf_iter 5000 -inbreed_iter 5000 -inbreed_tole 1e-9 \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/az_pop_mqu.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/az_pop_mqu \
&> /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/az_pop_mqu.log

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -inbreed 1 -inbreedSites -iter 5000 -maf_iter 5000 -inbreed_iter 5000 -inbreed_tole 1e-9 \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/tx_pop_mqu.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/tx_pop_mqu \
&> /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/tx_pop_mqu.log

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -inbreed 1 -inbreedSites -iter 5000 -maf_iter 5000 -inbreed_iter 5000 -inbreed_tole 1e-9 \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/nm_pop_mqu.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/nm_pop_mqu \
&> /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/nm_pop_mqu.log

#PCA using inbreeding

python /scratch/snyder/m/mathur20/MQU/2019/softwares/pcangsd/pcangsd.py -threads 80 \
-maf_tole 1e-9 -tole 1e-9 -iter 5000 -maf_iter 5000 \
-hwe /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/all_pop_mqu.lrt.sites.npy \
-beagle /scratch/snyder/m/mathur20/MQU/2020/angsd/gl/all_pop_mqu.beagle.gz \
-o /scratch/snyder/m/mathur20/MQU/2020/angsd/inbreeding/all_pop_mqu_pca