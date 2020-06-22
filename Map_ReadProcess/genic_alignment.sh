##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     genic_align.sh                        			###
###   To extract MQU genic reads and mapping to MQU referece genes		### 
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N align_genic


cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/3.6.0
module load samtools
module load BBMap/37.93

bwa index /scratch/snyder/m/mathur20/MQU/2020/adapt/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta
samtools faidx /scratch/snyder/m/mathur20/MQU/2020/adapt/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta

# Extract genic reads using bbmap #
while read -a line
do 
	bbsplit.sh -Xmx200g \
	minratio=0.85 minhits=2 \
	ref=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
	in=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/reads/${line[0]}_R1_filtered.fastq \
	in2=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/reads/${line[0]}_R2_filtered.fastq \
	out=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/reads/${line[0]}_genic.fastq \
	scafstats=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/stats/${line[0]}_bbstats.txt
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

#Aligning genic reads to MQU genes (Mathur et al. 2019) #

while read -a line
do 
	cd /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/bam/
	bwa mem -t 80 -M -R "@RG\tID:group1\tSM:${line[0]}\tPL:illumina\tLB:lib1\tPU:unit1" \
	/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
	/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/reads/${line[0]}_genic.fastq \
	> ${line[0]}_MQU_genic.sam
	PicardCommandLine ValidateSamFile I=${line[0]}_MQU_genic.sam MODE=SUMMARY O=${line[0]}_MQU_genic.sam.txt
	PicardCommandLine SortSam INPUT=${line[0]}_MQU_genic.sam OUTPUT=sorted_${line[0]}_MQU_genic.bam SORT_ORDER=coordinate
	PicardCommandLine MarkDuplicates INPUT=sorted_${line[0]}_MQU_genic.bam OUTPUT=dedup_${line[0]}_MQU_genic.bam METRICS_FILE=metrics_${line[0]}_MQU_genic.bam.txt
	PicardCommandLine BuildBamIndex INPUT=dedup_${line[0]}_MQU_genic.bam
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

PicardCommandLine CreateSequenceDictionary \
reference=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
output=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.dict

while read -a line
do
	GenomeAnalysisTK -nt 80 -T RealignerTargetCreator \
	-R /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
	-I /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/bam/dedup_${line[0]}_MQU_genic.bam \
	-o /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/gatk/forIndelRealigner.${line[0]}.intervals

	GenomeAnalysisTK -T IndelRealigner \
	-R /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/MQU_male_k60_SE-abyss_min5000.all.maker.transcripts.fasta \
	-I /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/bam/dedup_${line[0]}_MQU_genic.bam \
	-targetIntervals /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/gatk/forIndelRealigner.${line[0]}.intervals \
	-o /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/gatk/realigned_${line[0]}_reads.bam \
	&> /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/gatk/log/indelrealign.${line[0]}.logfile.txt
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

while read -a line
do
	PicardCommandLine FixMateInformation \
	INPUT=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/gatk/realigned_${line[0]}_reads.bam \
	OUTPUT=/scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/genic/gatk/final/${line[0]}.sorted.dedup.realigned.fixmate.bam \
	SO=coordinate \
	CREATE_INDEX=true
done < /scratch/snyder/m/mathur20/MQU/2019/try_final_angsd/align/sample.list

