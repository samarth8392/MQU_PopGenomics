
##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                    genetic_load.sh                       			###
###   			Estimating frequency of different impact SNPs			### 
###########################################################################


#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N genetic_load

cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load samtools
module load bedops
module load vcftools

# Extract SNPs from genes

awk '{gsub(/;$/,"");print}' /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/mqu_genes.gff \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/revisedgenome.mqu.gff

cat /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/revisedgenome.mqu.gff \
| awk '$3 =="gene" {print $0}' > /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/revisedgenome.mqu_justgenes.gff

gff2bed < /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/revisedgenome.mqu_justgenes.gff \
 > /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/MQU/revisedgenome.mqu_justgenes.bed 


vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/mqu/az/az_s2_eff.vcf --recode \
--stdout --bed /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/revisedgenome.all_justgenes.bed \
 > /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/r2/az_genic_snps_eff.vcf

vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/mqu/tx/tx_s2_eff.vcf --recode \
--stdout --bed /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/revisedgenome.all_justgenes.bed \
 > /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/r2/tx_genic_snps_eff.vcf

vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/snpeff/mqu/nm/nm_s2_eff.vcf --recode \
--stdout --bed /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/revisedgenome.all_justgenes.bed \
 > /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/r2/nm_genic_snps_eff.vcf

# Remove sites with warnings

grep "WARNING" az_s1_eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_az.txt
grep "WARNING" tx_s1_eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_tx.txt
grep "WARNING" nm_s1_eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_nm.txt


vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/az_genic_snps_eff.vcf --recode \
--stdout --exclude-positions /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/bad_sites_az.txt \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/az_genic_good_snps.vcf

vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/tx_genic_snps_eff.vcf --recode \
--stdout --exclude-positions /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/bad_sites_tx.txt \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/tx_genic_good_snps.vcf

vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/nm_genic_snps_eff.vcf --recode \
--stdout --exclude-positions /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/bad_sites_nm.txt \
> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/nm_genic_good_snps.vcf

# Get high, moderate, low, and no-impact variants

grep "HIGH" az_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > az_high_sites.txt
grep "HIGH" tx_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > tx_high_sites.txt
grep "HIGH" nm_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > nm_high_sites.txt

grep "LOW" az_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > az_low_sites.txt
grep "LOW" tx_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > tx_low_sites.txt
grep "LOW" nm_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > nm_low_sites.txt

grep "MODERATE" az_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > az_mod_sites.txt
grep "MODERATE" tx_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > tx_mod_sites.txt
grep "MODERATE" nm_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > nm_mod_sites.txt

grep "MODIFIER" az_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > az_nc_sites.txt
grep "MODIFIER" tx_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > tx_nc_sites.txt
grep "MODIFIER" nm_genic_good_snps.vcf | sed -e '1d' | cut -f1,2 > nm_nc_sites.txt

# Population Allele frequncies #
# Get frequency

for i in high low mod nc
do
	vcftools --bcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/az_s2.bcf --recode \
	--stdout --positions /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/r2/az_${i}_sites.txt \
	> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/r2/az_${i}_snps.vcf

	vcftools --bcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/tx_s2.bcf --recode \
	--stdout --positions /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/r2/tx_${i}_sites.txt \
	> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/r2/tx_${i}_snps.vcf

	vcftools --bcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/nm_s2.bcf --recode \
	--stdout --positions /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/r2/nm_${i}_sites.txt \
	> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/r2/nm_${i}_snps.vcf
done

for i in high low mod nc
do
	vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/r2/az_${i}_snps.vcf --freq \
	--out /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/r2/az_${i}_freq

	vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/r2/tx_${i}_snps.vcf --freq \
	--out /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/r2/tx_${i}_freq

	vcftools --vcf /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/r2/nm_${i}_snps.vcf --freq \
	--out /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/r2/nm_${i}_freq
done

# Get MAF

for i in high low mod nc
do
	sed -e '1d' /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/r2/az_${i}_freq.frq | cut -f6 \
	> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/az/r2/az_${i}_minor.freq
	sed -e '1d' /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/r2/tx_${i}_freq.frq | cut -f6 \
	> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/tx/r2/tx_${i}_minor.freq
	sed -e '1d' /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/r2/nm_${i}_freq.frq | cut -f6 \
	> /scratch/snyder/m/mathur20/MQU/2020/angsd/snps/heter/mqu/nm/r2/nm_${i}_minor.freq
done

## FOR ESTIMATING GENETIC LOAD FOR EACH INDIIVUDAL ##
# To get the allele counts from the genotypes #
# 0/0 = 0 allele ; 0/1 = 1 allele; 1/1/ = 2 alleles
for i in 1 2 3 4 5 6 7
do
	cut -f${i} az_genic_high.geno | cut -f1 -d ':' > ind_geno/az_ind${i}_high.geno
	cut -f${i} tx_genic_high.geno | cut -f1 -d ':' > ind_geno/tx_ind${i}_high.geno
	cut -f${i} nm_genic_high.geno | cut -f1 -d ':' > ind_geno/nm_ind${i}_high.geno

	cut -f${i} az_genic_moderate.geno | cut -f1 -d ':' > ind_geno/az_ind${i}_moderate.geno
	cut -f${i} tx_genic_moderate.geno | cut -f1 -d ':' > ind_geno/tx_ind${i}_moderate.geno
	cut -f${i} nm_genic_moderate.geno | cut -f1 -d ':' > ind_geno/nm_ind${i}_moderate.geno

	cut -f${i} az_genic_low.geno | cut -f1 -d ':' > ind_geno/az_ind${i}_low.geno
	cut -f${i} tx_genic_low.geno | cut -f1 -d ':' > ind_geno/tx_ind${i}_low.geno
	cut -f${i} nm_genic_low.geno | cut -f1 -d ':' > ind_geno/nm_ind${i}_low.geno

	cut -f${i} az_genic_modifier.geno | cut -f1 -d ':' > ind_geno/az_ind${i}_modifier.geno
	cut -f${i} tx_genic_modifier.geno | cut -f1 -d ':' > ind_geno/tx_ind${i}_modifier.geno
	cut -f${i} nm_genic_modifier.geno | cut -f1 -d ':' > ind_geno/nm_ind${i}_modifier.geno
done

pr -mts' ' az_*_high.geno > az_high_allind_allele
pr -mts' ' tx_*_high.geno > az_high_allind_allele
pr -mts' ' nm_*_high.geno > az_high_allind_allele
