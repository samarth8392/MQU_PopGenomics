##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     mu_rate.sh                       				###
###   		Estimating mutation rate in MQU								### 
###########################################################################

#!/bin/sh -l
#PBS -q debug
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=00:30:00
#PBS -N mu_rate

cd $PBS_O_WORKDIR
module load bioinfo
module load LASTZ/1.04.00
module load samtools
module load R/3.4.2

# MQU genome was aligned to chicken using LASTZ with 
# the conditions as follows: 
# T = 2, C = 2, H = 2,000, Y = 3,400, L = 6,000, K = 2,200. 

# Polymorphic loci were identified according to the following standards.
# 1) The nucleotide from either target or query was not classified as a N or n;
# 2) The locus was not in an alignment gap; 

# The final mutation rate per nt per year µwas calculated with the following formula:
# µ= (counts of mutated loci / sequence length) / 2t. The divergence time t

lastz_32 /scratch/snyder/m/mathur20/MQU/2019/angsd/final_74/ref_chick/bam/chicken_genome.fa[multiple] \
/scratch/snyder/m/mathur20/MQU/2019/angsd/final_74/ref_MQU/folded/MQU_male.min500.fa \
--ydrop=3400 --notransition --gapped --step=20 \
--inner=2000 --gappedthresh=6000 --hspthresh=2200 \
--allocate:traceback=10485760 --allocate:target=1048576000 --allocate:query=1048576000 \
--format=general:name1,size1,name2,size2,nmatch,nmismatch,ngap,cgap \
--ambiguous=iupac > /scratch/snyder/m/mathur20/MQU/2019/msmc/lastz/lastz_chick_mqu.tab

samtools view -bS lastz_chick_mqu.sam > lastz_chick_mqu.bam

samtools index lastz_chick_mqu.bam 


# Get the number of mismatches
Rscript /scratch/snyder/m/mathur20/MQU/2019/msmc/jobcodes/mu_rate.R \
-f /scratch/snyder/m/mathur20/MQU/2019/msmc/lastz/lastz_chick_mqu.tab \
-o /scratch/snyder/m/mathur20/MQU/2019/msmc/lastz/mu_rate_allchr.txt