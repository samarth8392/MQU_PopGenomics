###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     mito_haplotypes.sh                        		### 
### To get SNPs, filter SNPs, and get concensus mito haplotypes from    ###
### reads mapped to the trimmed reference sequence              		###
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N mito_haps

#cd $PBS_O_WORKDIR
#module purge


# STEPS 1-7 : See mitogenome_assembly.sh #

# STEP 8: Get SNPs *

cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_jobs/
for prefix in $(ls *.job | sed -r 's/_step[012345678][.]job//' | uniq)
do
 echo  "#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N ${prefix}_Step8
cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load samtools
module load bcftools

samtools faidx /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female_mito_nonumts_trim2.fasta

samtools mpileup -t DP -t SP -u -g -C 75 \
-f /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female_mito_nonumts_trim2.fasta \
/scratch/snyder/m/mathur20/MQU/2019/mito/round3/alignments/dedup_sorted_${prefix}_MQU_nonumt_mito.bam \
| bcftools call -c -v > /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_haps/${prefix}_female_mito.snps.raw.vcf"> ${prefix}_step8.job
done

# STEP 9: Filter SNPs and get haplotypes #

cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_jobs/
for prefix in $(ls *.job | sed -r 's/_step[0123456789][.]job//' | uniq)
do
 echo  "#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N ${prefix}_Step9
cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load samtools
module load bcftools
module load vcflib

#cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_haps

vcffilter -f \"DP > 10\" \
/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_haps/${prefix}_female_mito.snps.raw.vcf > ${prefix}_female_mito.snps.filter.vcf

bgzip /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_haps/${prefix}_female_mito.snps.filter.vcf

bcftools index /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_haps/${prefix}_female_mito.snps.filter.vcf.gz

cat /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female_mito_nonumts_trim2.fasta \
| bcftools consensus /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_haps/${prefix}_female_mito.snps.filter.vcf.gz \
> /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_haps/${prefix}_mito_haplotype.fa"> ${prefix}_step9.job
done

# STEP 9b: Rename haplotypes #

cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/fasta/
for prefix in $(ls *.fa | sed -r 's/_mito[.]fa//' | uniq)
do
 cat ${prefix}_mito.fa | awk -v a="$prefix" '/^>/{print ">" a " Mitochondrion, complete genome "; next}{print}' > ${prefix}_mito_renamed.fa
done

################## RUN JOB CODES ############################


cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_jobs/
for prefix in $(ls *.job | sed -r 's/_step[0123456789][.]job//' | uniq)
do
 qsub ${prefix}_step9.job
done
