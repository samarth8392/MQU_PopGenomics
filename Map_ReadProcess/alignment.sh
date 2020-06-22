##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     align.sh                        				###
###   To align the MQU WGR reads to the reference assemblies			### 
###########################################################################

#!/bin/sh -l
#PBS -q fnrquail
#PBS -l nodes=1:ppn=5,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N align_depth_final


cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/3.6.0
module load samtools

# get reads
while read -a line
do 
	cp /scratch/snyder/m/mathur20/MQU/2019/reads/final/${line[0]}* \
	/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/bam/reads/
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/bam/sample.list


# Align to MQU male assembly (Can change it to chicken for chicken BAM)

bwa index /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa
samtools faidx /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa

while read -a line
do 
	cd /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam/
	bwa mem -t 80 -M -R "@RG\tID:group1\tSM:${line[0]}\tPL:illumina\tLB:lib1\tPU:unit1" \
	/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
	/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/reads/${line[0]}_R1_filtered.fastq \
	/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/reads/${line[0]}_R2_filtered.fastq \
	> ${line[0]}_MQU_male.sam
	PicardCommandLine ValidateSamFile I=${line[0]}_MQU_male.sam MODE=SUMMARY O=${line[0]}_MQU_male.sam.txt
	PicardCommandLine SortSam INPUT=${line[0]}_MQU_male.sam OUTPUT=sorted_${line[0]}_MQU_male.bam SORT_ORDER=coordinate
	PicardCommandLine MarkDuplicates INPUT=sorted_${line[0]}_MQU_male.bam OUTPUT=dedup_${line[0]}_MQU_male.bam METRICS_FILE=metrics_${line[0]}_MQU_male.bam.txt
	PicardCommandLine BuildBamIndex INPUT=dedup_${line[0]}_MQU_male.bam
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

# Get coverage depth and breadth
while read -a line
do
	samtools depth -a /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam/dedup_${line[0]}_MQU_male.bam \
	| awk '{c++;s+=$3}END{print s/c}' \
	> /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam_stats/${line[0]}_stats.txt

	samtools depth -a /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam/dedup_${line[0]}_MQU_male.bam \
	| awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' \
	> /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam_stats/${line[0]}_1x_breadth.txt
	
	samtools flagstat /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam/dedup_${line[0]}_MQU_male.bam \
	| awk -F "[(|%]" 'NR == 3 {print $2}' \
	> /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam_stats/${line[0]}_mapped.txt
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

# Realign indels

PicardCommandLine CreateSequenceDictionary \
reference=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
output=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.dict

while read -a line
do
	GenomeAnalysisTK -nt 80 -T RealignerTargetCreator \
	-R /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
	-I /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam/dedup_${line[0]}_MQU_male.bam \
	-o /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/forIndelRealigner.${line[0]}.intervals

	GenomeAnalysisTK -T IndelRealigner \
	-R /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
	-I /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam/dedup_${line[0]}_MQU_male.bam \
	-targetIntervals /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/forIndelRealigner.${line[0]}.intervals \
	-o /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/realigned_${line[0]}_reads.bam \
	&> /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/log/indelrealign.${line[0]}.logfile.txt
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

# Fix mate pair info in BAM
while read -a line
do
	PicardCommandLine FixMateInformation \
	INPUT=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/realigned_${line[0]}_reads.bam \
	OUTPUT=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/${line[0]}.sorted.dedup.realigned.fixmate.bam \
	SO=coordinate \
	CREATE_INDEX=true
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

#Base quality score recalibration
while read -a line
do
	GenomeAnalysisTK -T BaseRecalibrator \
	-R /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
	-knownSites /scratch/snyder/m/mathur20/MQU/2018/markers/neutral/snps/raw_indels.vcf \
	-knownSites /scratch/snyder/m/mathur20/MQU/2018/markers/neutral/snps/raw_snps.vcf \
	-I /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/${line[0]}.sorted.dedup.realigned.fixmate.bam \
	-o /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/table/${line[0]}.sorted.dedup.realigned.fixmate.recal_data.table \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate

	GenomeAnalysisTK -T PrintReads \
	-R /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/MQU_male.min500.fa \
	-BQSR /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/table/${line[0]}.sorted.dedup.realigned.fixmate.recal_data.table \
	-I /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/${line[0]}.sorted.dedup.realigned.fixmate.bam \
	-o /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/final/${line[0]}.sorted.dedup.realigned.fixmate.recal.bam
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

# Get final coverage depth and breadth
while read -a line
do
	samtools depth -a /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/final/${line[0]}.sorted.dedup.realigned.fixmate.recal.bam \
	| awk '{c++;s+=$3}END{print s/c}' \
	> /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam_stats/final/${line[0]}_stats.txt

	samtools depth -a /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/gatk/final/${line[0]}.sorted.dedup.realigned.fixmate.recal.bam \
	| awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' \
	> /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam_stats/final/${line[0]}_1x_breadth.txt
	
	samtools flagstat /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam/dedup_${line[0]}_MQU_male.bam \
	> /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/bam_stats/final/${line[0]}_mapped.txt
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list
