###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     SNPeff.sh      			                  		###
###   To functionally annotate SNP markers 						    	###
###                                                                     ###
###########################################################################

#!/bin/sh -l
#PBS -q debug
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=00:30:00
#PBS -N SNPeff

cd $PBS_O_WORKDIR

module use /apps/group/bioinformatics/modules
module load snpEff/4.3
module load bedops

#snpEff build -c /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/snpEff.config -gff3 -v MQU &> build.logfile.txt

snpEff ann -stats -c /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/snpEff.config \
-no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v \
MQU /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/az_s1.vcf \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/az/az_s1_eff.vcf

snpEff ann -stats -c /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/snpEff.config \
-no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v \
MQU /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/tx_s1.vcf \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/tx/tx_s1_eff.vcf

snpEff ann -stats -c /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/snpEff.config \
-no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v \
MQU /scratch/snyder/m/mathur20/MQU/2020/angsd/subsample/s1/nm_s1.vcf \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/nm/nm_s1_eff.vcf

# vcf to bed
vcf2bed < /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/az/az_s1_eff.vcf > /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/az/SNPeff.az_mqu_wgs.bed 
vcf2bed < /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/tx/tx_s1_eff.vcf > /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/tx/SNPeff.tx_mqu_wgs.bed 
vcf2bed < /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/nm/nm_s1_eff.vcf > /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/nm/SNPeff.nm_mqu_wgs.bed 