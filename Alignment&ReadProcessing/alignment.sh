##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 03/29/19                  Last Modified: 03/29/19 ###
###########################################################################
###########################################################################
###                     alignment.sh                        			###
###   To align the MQU WGR reads to the reference assemblies			### 
###########################################################################

#!/bin/sh -l
#PBS -q fnrdewoody
#PBS -l nodes=1:ppn=20,naccesspolicy=shared
#PBS -l walltime=336:00:00
#PBS -N align_reads


cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/3.6.0
module load samtools


##### FIRST WAY: ALL TOGETHER #####

# index the fasta files

#bwa index MQU_male.min500.fa
#samtools faidx MQU_male.min500.fa

#bwa index MQU_female.min500.fasta
#samtools faidx MQU_female.min500.fasta

# Concatenate the fastq files

cat /scratch/snyder/m/mathur20/MQU/2019/reads/run2/*_R1_filtered.fastq > /scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/all_wgr_R1.fastq
cat /scratch/snyder/m/mathur20/MQU/2019/reads/run2/*_R2_filtered.fastq > /scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/all_wgr_R2.fastq

# map reads to the reference using BWA

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:all_wgr\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/all_wgr_R1.fastq \
/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/all_wgr_R2.fastq > /scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/MQU_male_all_wgr.sam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_female.min500.fasta \
/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/all_wgr_R1.fastq \
/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/all_wgr_R2.fastq > /scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/MQU_female_all_wgr.sam

# Validate the SAM file should produce a validate_output.txt file that says there are no errors.

PicardCommandLine ValidateSamFile I=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/MQU_male_all_wgr.sam MODE=SUMMARY O=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/MQU_male_all_wgr.sam.txt
PicardCommandLine ValidateSamFile I=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/MQU_female_all_wgr.sam MODE=SUMMARY O=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/MQU_female_all_wgr.sam_samfile.txt

# Create BAM files

PicardCommandLine SortSam INPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/MQU_male_all_wgr.sam OUTPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/sorted_MQU_male_all_wgr.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/sorted_MQU_male_all_wgr.bam OUTPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/dedup_MQU_male_all_wgr.bam METRICS_FILE=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/metrics_MQU_male_all_wgr.txt
PicardCommandLine BuildBamIndex INPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/dedup_MQU_male_all_wgr.bam

PicardCommandLine SortSam INPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/MQU_female_all_wgr.sam OUTPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/sorted_MQU_female_all_wgr.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/sorted_MQU_female_all_wgr.bam OUTPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/dedup_MQU_female_all_wgr.bam METRICS_FILE=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/metrics_MQU_female_all_wgr.txt
PicardCommandLine BuildBamIndex INPUT=/scratch/snyder/m/mathur20/MQU/2019/alignments/all_wgr/dedup_MQU_female_all_wgr.bam

##### END OF FIRST WAY: ALL TOGETHER #####

##### SECOND WAY: INDIVIDUALLY #####

# Aligning to the male assembly

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6537\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037268_E6537_S32_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037268_E6537_S32_R2_filtered.fastq > MQU_male_E6537.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6537.sam MODE=SUMMARY O=MQU_male_E6537.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6537.sam OUTPUT=sorted_MQU_male_E6537.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6537.bam OUTPUT=dedup_MQU_male_E6537.bam METRICS_FILE=metrics_MQU_male_E6537.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6537.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6538\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037269_E6538_S33_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037269_E6538_S33_R2_filtered.fastq > MQU_male_E6538.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6538.sam MODE=SUMMARY O=MQU_male_E6538.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6538.sam OUTPUT=sorted_MQU_male_E6538.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6538.bam OUTPUT=dedup_MQU_male_E6538.bam METRICS_FILE=metrics_MQU_male_E6538.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6538.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6539\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037270_E6539_S34_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037270_E6539_S34_R2_filtered.fastq > MQU_male_E6539.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6539.sam MODE=SUMMARY O=MQU_male_E6539.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6539.sam OUTPUT=sorted_MQU_male_E6539.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6539.bam OUTPUT=dedup_MQU_male_E6539.bam METRICS_FILE=metrics_MQU_male_E6539.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6539.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6541\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037271_E6541_S35_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037271_E6541_S35_R2_filtered.fastq > MQU_male_E6541.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6541.sam MODE=SUMMARY O=MQU_male_E6541.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6541.sam OUTPUT=sorted_MQU_male_E6541.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6541.bam OUTPUT=dedup_MQU_male_E6541.bam METRICS_FILE=metrics_MQU_male_E6541.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6541.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6542\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037272_E6542_S5_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037272_E6542_S5_R2_filtered.fastq > MQU_male_E6542.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6542.sam MODE=SUMMARY O=MQU_male_E6542.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6542.sam OUTPUT=sorted_MQU_male_E6542.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6542.bam OUTPUT=dedup_MQU_male_E6542.bam METRICS_FILE=metrics_MQU_male_E6542.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6542.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6543\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037273_E6543_S37_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037273_E6543_S37_R2_filtered.fastq > MQU_male_E6543.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6543.sam MODE=SUMMARY O=MQU_male_E6543.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6543.sam OUTPUT=sorted_MQU_male_E6543.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6543.bam OUTPUT=dedup_MQU_male_E6543.bam METRICS_FILE=metrics_MQU_male_E6543.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6543.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6545\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037274_E6545_S38_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037274_E6545_S38_R2_filtered.fastq > MQU_male_E6545.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6545.sam MODE=SUMMARY O=MQU_male_E6545.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6545.sam OUTPUT=sorted_MQU_male_E6545.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6545.bam OUTPUT=dedup_MQU_male_E6545.bam METRICS_FILE=metrics_MQU_male_E6545.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6545.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6548\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037275_E6548_S39_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037275_E6548_S39_R2_filtered.fastq > MQU_male_E6548.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6548.sam MODE=SUMMARY O=MQU_male_E6548.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6548.sam OUTPUT=sorted_MQU_male_E6548.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6548.bam OUTPUT=dedup_MQU_male_E6548.bam METRICS_FILE=metrics_MQU_male_E6548.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6548.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6566\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037276_E6566_S40_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037276_E6566_S40_R2_filtered.fastq > MQU_male_E6566.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6566.sam MODE=SUMMARY O=MQU_male_E6566.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6566.sam OUTPUT=sorted_MQU_male_E6566.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6566.bam OUTPUT=dedup_MQU_male_E6566.bam METRICS_FILE=metrics_MQU_male_E6566.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6566.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6597\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037277_E6597_S41_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037277_E6597_S41_R2_filtered.fastq > MQU_male_E6597.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6597.sam MODE=SUMMARY O=MQU_male_E6597.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6597.sam OUTPUT=sorted_MQU_male_E6597.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6597.bam OUTPUT=dedup_MQU_male_E6597.bam METRICS_FILE=metrics_MQU_male_E6597.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6597.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6598\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037278_E6598_S42_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037278_E6598_S42_R2_filtered.fastq > MQU_male_E6598.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6598.sam MODE=SUMMARY O=MQU_male_E6598.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6598.sam OUTPUT=sorted_MQU_male_E6598.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6598.bam OUTPUT=dedup_MQU_male_E6598.bam METRICS_FILE=metrics_MQU_male_E6598.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6598.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6599\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037279_E6599_S43_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037279_E6599_S43_R2_filtered.fastq > MQU_male_E6599.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6599.sam MODE=SUMMARY O=MQU_male_E6599.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6599.sam OUTPUT=sorted_MQU_male_E6599.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6599.bam OUTPUT=dedup_MQU_male_E6599.bam METRICS_FILE=metrics_MQU_male_E6599.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6599.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6600\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037280_E6600_S44_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037280_E6600_S44_R2_filtered.fastq > MQU_male_E6600.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6600.sam MODE=SUMMARY O=MQU_male_E6600.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6600.sam OUTPUT=sorted_MQU_male_E6600.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6600.bam OUTPUT=dedup_MQU_male_E6600.bam METRICS_FILE=metrics_MQU_male_E6600.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6600.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6609\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037281_E6609_S45_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037281_E6609_S45_R2_filtered.fastq > MQU_male_E6609.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6609.sam MODE=SUMMARY O=MQU_male_E6609.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6609.sam OUTPUT=sorted_MQU_male_E6609.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6609.bam OUTPUT=dedup_MQU_male_E6609.bam METRICS_FILE=metrics_MQU_male_E6609.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6609.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6612\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037282_E6612_S6_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037282_E6612_S6_R2_filtered.fastq > MQU_male_E6612.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6612.sam MODE=SUMMARY O=MQU_male_E6612.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6612.sam OUTPUT=sorted_MQU_male_E6612.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6612.bam OUTPUT=dedup_MQU_male_E6612.bam METRICS_FILE=metrics_MQU_male_E6612.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6612.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6705\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037283_E6705_S47_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037283_E6705_S47_R2_filtered.fastq > MQU_male_E6705.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6705.sam MODE=SUMMARY O=MQU_male_E6705.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6705.sam OUTPUT=sorted_MQU_male_E6705.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6705.bam OUTPUT=dedup_MQU_male_E6705.bam METRICS_FILE=metrics_MQU_male_E6705.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6705.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6715\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037284_E6715_S48_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037284_E6715_S48_R2_filtered.fastq > MQU_male_E6715.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6715.sam MODE=SUMMARY O=MQU_male_E6715.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6715.sam OUTPUT=sorted_MQU_male_E6715.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6715.bam OUTPUT=dedup_MQU_male_E6715.bam METRICS_FILE=metrics_MQU_male_E6715.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6715.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6716\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037285_E6716_S49_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037285_E6716_S49_R2_filtered.fastq > MQU_male_E6716.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6716.sam MODE=SUMMARY O=MQU_male_E6716.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6716.sam OUTPUT=sorted_MQU_male_E6716.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6716.bam OUTPUT=dedup_MQU_male_E6716.bam METRICS_FILE=metrics_MQU_male_E6716.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6716.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6718\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037286_E6718_S50_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037286_E6718_S50_R2_filtered.fastq > MQU_male_E6718.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6718.sam MODE=SUMMARY O=MQU_male_E6718.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6718.sam OUTPUT=sorted_MQU_male_E6718.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6718.bam OUTPUT=dedup_MQU_male_E6718.bam METRICS_FILE=metrics_MQU_male_E6718.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6718.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6758\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037287_E6758_S51_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037287_E6758_S51_R2_filtered.fastq > MQU_male_E6758.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6758.sam MODE=SUMMARY O=MQU_male_E6758.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6758.sam OUTPUT=sorted_MQU_male_E6758.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6758.bam OUTPUT=dedup_MQU_male_E6758.bam METRICS_FILE=metrics_MQU_male_E6758.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6758.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6759\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037288_E6759_S52_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037288_E6759_S52_R2_filtered.fastq > MQU_male_E6759.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6759.sam MODE=SUMMARY O=MQU_male_E6759.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6759.sam OUTPUT=sorted_MQU_male_E6759.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6759.bam OUTPUT=dedup_MQU_male_E6759.bam METRICS_FILE=metrics_MQU_male_E6759.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6759.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6762\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037289_E6762_S53_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037289_E6762_S53_R2_filtered.fastq > MQU_male_E6762.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6762.sam MODE=SUMMARY O=MQU_male_E6762.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6762.sam OUTPUT=sorted_MQU_male_E6762.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6762.bam OUTPUT=dedup_MQU_male_E6762.bam METRICS_FILE=metrics_MQU_male_E6762.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6762.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6799\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037290_E6799_S54_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037290_E6799_S54_R2_filtered.fastq > MQU_male_E6799.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6799.sam MODE=SUMMARY O=MQU_male_E6799.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6799.sam OUTPUT=sorted_MQU_male_E6799.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6799.bam OUTPUT=dedup_MQU_male_E6799.bam METRICS_FILE=metrics_MQU_male_E6799.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6799.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6800\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037291_E6800_S55_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037291_E6800_S55_R2_filtered.fastq > MQU_male_E6800.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6800.sam MODE=SUMMARY O=MQU_male_E6800.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6800.sam OUTPUT=sorted_MQU_male_E6800.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6800.bam OUTPUT=dedup_MQU_male_E6800.bam METRICS_FILE=metrics_MQU_male_E6800.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6800.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6801\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037292_E6801_S56_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037292_E6801_S56_R2_filtered.fastq > MQU_male_E6801.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6801.sam MODE=SUMMARY O=MQU_male_E6801.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6801.sam OUTPUT=sorted_MQU_male_E6801.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6801.bam OUTPUT=dedup_MQU_male_E6801.bam METRICS_FILE=metrics_MQU_male_E6801.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6801.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6803\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037293_E6803_S57_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037293_E6803_S57_R2_filtered.fastq > MQU_male_E6803.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6803.sam MODE=SUMMARY O=MQU_male_E6803.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6803.sam OUTPUT=sorted_MQU_male_E6803.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6803.bam OUTPUT=dedup_MQU_male_E6803.bam METRICS_FILE=metrics_MQU_male_E6803.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6803.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6841\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037294_E6841_S100_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037294_E6841_S100_R2_filtered.fastq > MQU_male_E6841.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6841.sam MODE=SUMMARY O=MQU_male_E6841.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6841.sam OUTPUT=sorted_MQU_male_E6841.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6841.bam OUTPUT=dedup_MQU_male_E6841.bam METRICS_FILE=metrics_MQU_male_E6841.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6841.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6843\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037295_E6843_S101_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037295_E6843_S101_R2_filtered.fastq > MQU_male_E6843.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6843.sam MODE=SUMMARY O=MQU_male_E6843.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6843.sam OUTPUT=sorted_MQU_male_E6843.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6843.bam OUTPUT=dedup_MQU_male_E6843.bam METRICS_FILE=metrics_MQU_male_E6843.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6843.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6847\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037296_E6847_S60_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037296_E6847_S60_R2_filtered.fastq > MQU_male_E6847.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6847.sam MODE=SUMMARY O=MQU_male_E6847.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6847.sam OUTPUT=sorted_MQU_male_E6847.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6847.bam OUTPUT=dedup_MQU_male_E6847.bam METRICS_FILE=metrics_MQU_male_E6847.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6847.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6960\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037297_E6960_S102_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037297_E6960_S102_R2_filtered.fastq > MQU_male_E6960.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6960.sam MODE=SUMMARY O=MQU_male_E6960.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6960.sam OUTPUT=sorted_MQU_male_E6960.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6960.bam OUTPUT=dedup_MQU_male_E6960.bam METRICS_FILE=metrics_MQU_male_E6960.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6960.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6963\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037298_E6963_S103_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037298_E6963_S103_R2_filtered.fastq > MQU_male_E6963.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6963.sam MODE=SUMMARY O=MQU_male_E6963.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6963.sam OUTPUT=sorted_MQU_male_E6963.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6963.bam OUTPUT=dedup_MQU_male_E6963.bam METRICS_FILE=metrics_MQU_male_E6963.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6963.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6964\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037299_E6964_S63_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037299_E6964_S63_R2_filtered.fastq > MQU_male_E6964.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6964.sam MODE=SUMMARY O=MQU_male_E6964.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6964.sam OUTPUT=sorted_MQU_male_E6964.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6964.bam OUTPUT=dedup_MQU_male_E6964.bam METRICS_FILE=metrics_MQU_male_E6964.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6964.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E6966\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037300_E6966_S64_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037300_E6966_S64_R2_filtered.fastq > MQU_male_E6966.sam
PicardCommandLine ValidateSamFile I=MQU_male_E6966.sam MODE=SUMMARY O=MQU_male_E6966.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E6966.sam OUTPUT=sorted_MQU_male_E6966.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E6966.bam OUTPUT=dedup_MQU_male_E6966.bam METRICS_FILE=metrics_MQU_male_E6966.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E6966.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7010\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037301_E7010_S104_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037301_E7010_S104_R2_filtered.fastq > MQU_male_E7010.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7010.sam MODE=SUMMARY O=MQU_male_E7010.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7010.sam OUTPUT=sorted_MQU_male_E7010.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7010.bam OUTPUT=dedup_MQU_male_E7010.bam METRICS_FILE=metrics_MQU_male_E7010.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7010.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7011\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037302_E7011_S105_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037302_E7011_S105_R2_filtered.fastq > MQU_male_E7011.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7011.sam MODE=SUMMARY O=MQU_male_E7011.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7011.sam OUTPUT=sorted_MQU_male_E7011.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7011.bam OUTPUT=dedup_MQU_male_E7011.bam METRICS_FILE=metrics_MQU_male_E7011.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7011.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7012\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037303_E7012_S106_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037303_E7012_S106_R2_filtered.fastq > MQU_male_E7012.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7012.sam MODE=SUMMARY O=MQU_male_E7012.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7012.sam OUTPUT=sorted_MQU_male_E7012.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7012.bam OUTPUT=dedup_MQU_male_E7012.bam METRICS_FILE=metrics_MQU_male_E7012.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7012.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7035\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037304_E7035_S68_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037304_E7035_S68_R2_filtered.fastq > MQU_male_E7035.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7035.sam MODE=SUMMARY O=MQU_male_E7035.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7035.sam OUTPUT=sorted_MQU_male_E7035.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7035.bam OUTPUT=dedup_MQU_male_E7035.bam METRICS_FILE=metrics_MQU_male_E7035.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7035.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7036\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037305_E7036_S107_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037305_E7036_S107_R2_filtered.fastq > MQU_male_E7036.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7036.sam MODE=SUMMARY O=MQU_male_E7036.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7036.sam OUTPUT=sorted_MQU_male_E7036.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7036.bam OUTPUT=dedup_MQU_male_E7036.bam METRICS_FILE=metrics_MQU_male_E7036.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7036.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7048\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037306_E7048_S108_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037306_E7048_S108_R2_filtered.fastq > MQU_male_E7048.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7048.sam MODE=SUMMARY O=MQU_male_E7048.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7048.sam OUTPUT=sorted_MQU_male_E7048.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7048.bam OUTPUT=dedup_MQU_male_E7048.bam METRICS_FILE=metrics_MQU_male_E7048.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7048.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7049\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037307_E7049_S109_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037307_E7049_S109_R2_filtered.fastq > MQU_male_E7049.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7049.sam MODE=SUMMARY O=MQU_male_E7049.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7049.sam OUTPUT=sorted_MQU_male_E7049.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7049.bam OUTPUT=dedup_MQU_male_E7049.bam METRICS_FILE=metrics_MQU_male_E7049.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7049.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7159\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037308_E7159_S110_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037308_E7159_S110_R2_filtered.fastq > MQU_male_E7159.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7159.sam MODE=SUMMARY O=MQU_male_E7159.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7159.sam OUTPUT=sorted_MQU_male_E7159.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7159.bam OUTPUT=dedup_MQU_male_E7159.bam METRICS_FILE=metrics_MQU_male_E7159.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7159.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7160\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037309_E7160_S111_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037309_E7160_S111_R2_filtered.fastq > MQU_male_E7160.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7160.sam MODE=SUMMARY O=MQU_male_E7160.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7160.sam OUTPUT=sorted_MQU_male_E7160.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7160.bam OUTPUT=dedup_MQU_male_E7160.bam METRICS_FILE=metrics_MQU_male_E7160.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7160.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7161\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037310_E7161_S112_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037310_E7161_S112_R2_filtered.fastq > MQU_male_E7161.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7161.sam MODE=SUMMARY O=MQU_male_E7161.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7161.sam OUTPUT=sorted_MQU_male_E7161.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7161.bam OUTPUT=dedup_MQU_male_E7161.bam METRICS_FILE=metrics_MQU_male_E7161.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7161.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7175\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037311_E7175_S113_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037311_E7175_S113_R2_filtered.fastq > MQU_male_E7175.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7175.sam MODE=SUMMARY O=MQU_male_E7175.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7175.sam OUTPUT=sorted_MQU_male_E7175.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7175.bam OUTPUT=dedup_MQU_male_E7175.bam METRICS_FILE=metrics_MQU_male_E7175.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7175.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7177\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037312_E7177_S114_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037312_E7177_S114_R2_filtered.fastq > MQU_male_E7177.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7177.sam MODE=SUMMARY O=MQU_male_E7177.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7177.sam OUTPUT=sorted_MQU_male_E7177.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7177.bam OUTPUT=dedup_MQU_male_E7177.bam METRICS_FILE=metrics_MQU_male_E7177.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7177.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7219\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037313_E7219_S77_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037313_E7219_S77_R2_filtered.fastq > MQU_male_E7219.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7219.sam MODE=SUMMARY O=MQU_male_E7219.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7219.sam OUTPUT=sorted_MQU_male_E7219.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7219.bam OUTPUT=dedup_MQU_male_E7219.bam METRICS_FILE=metrics_MQU_male_E7219.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7219.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7755\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037314_E7755_S78_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037314_E7755_S78_R2_filtered.fastq > MQU_male_E7755.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7755.sam MODE=SUMMARY O=MQU_male_E7755.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7755.sam OUTPUT=sorted_MQU_male_E7755.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7755.bam OUTPUT=dedup_MQU_male_E7755.bam METRICS_FILE=metrics_MQU_male_E7755.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7755.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7757\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037315_E7757_S79_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037315_E7757_S79_R2_filtered.fastq > MQU_male_E7757.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7757.sam MODE=SUMMARY O=MQU_male_E7757.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7757.sam OUTPUT=sorted_MQU_male_E7757.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7757.bam OUTPUT=dedup_MQU_male_E7757.bam METRICS_FILE=metrics_MQU_male_E7757.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7757.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7879\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037316_E7879_S115_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037316_E7879_S115_R2_filtered.fastq > MQU_male_E7879.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7879.sam MODE=SUMMARY O=MQU_male_E7879.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7879.sam OUTPUT=sorted_MQU_male_E7879.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7879.bam OUTPUT=dedup_MQU_male_E7879.bam METRICS_FILE=metrics_MQU_male_E7879.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7879.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7910\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037317_E7910_S116_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037317_E7910_S116_R2_filtered.fastq > MQU_male_E7910.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7910.sam MODE=SUMMARY O=MQU_male_E7910.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7910.sam OUTPUT=sorted_MQU_male_E7910.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7910.bam OUTPUT=dedup_MQU_male_E7910.bam METRICS_FILE=metrics_MQU_male_E7910.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7910.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7915\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037318_E7915_S117_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037318_E7915_S117_R2_filtered.fastq > MQU_male_E7915.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7915.sam MODE=SUMMARY O=MQU_male_E7915.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7915.sam OUTPUT=sorted_MQU_male_E7915.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7915.bam OUTPUT=dedup_MQU_male_E7915.bam METRICS_FILE=metrics_MQU_male_E7915.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7915.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E7940\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037319_E7940_S118_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037319_E7940_S118_R2_filtered.fastq > MQU_male_E7940.sam
PicardCommandLine ValidateSamFile I=MQU_male_E7940.sam MODE=SUMMARY O=MQU_male_E7940.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E7940.sam OUTPUT=sorted_MQU_male_E7940.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E7940.bam OUTPUT=dedup_MQU_male_E7940.bam METRICS_FILE=metrics_MQU_male_E7940.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E7940.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8002\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037320_E8002_S84_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037320_E8002_S84_R2_filtered.fastq > MQU_male_E8002.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8002.sam MODE=SUMMARY O=MQU_male_E8002.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8002.sam OUTPUT=sorted_MQU_male_E8002.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8002.bam OUTPUT=dedup_MQU_male_E8002.bam METRICS_FILE=metrics_MQU_male_E8002.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8002.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8035\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037321_E8035_S100_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037321_E8035_S100_R2_filtered.fastq > MQU_male_E8035.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8035.sam MODE=SUMMARY O=MQU_male_E8035.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8035.sam OUTPUT=sorted_MQU_male_E8035.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8035.bam OUTPUT=dedup_MQU_male_E8035.bam METRICS_FILE=metrics_MQU_male_E8035.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8035.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8036\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037322_E8036_S86_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037322_E8036_S86_R2_filtered.fastq > MQU_male_E8036.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8036.sam MODE=SUMMARY O=MQU_male_E8036.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8036.sam OUTPUT=sorted_MQU_male_E8036.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8036.bam OUTPUT=dedup_MQU_male_E8036.bam METRICS_FILE=metrics_MQU_male_E8036.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8036.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8063\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037323_E8063_S101_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037323_E8063_S101_R2_filtered.fastq > MQU_male_E8063.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8063.sam MODE=SUMMARY O=MQU_male_E8063.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8063.sam OUTPUT=sorted_MQU_male_E8063.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8063.bam OUTPUT=dedup_MQU_male_E8063.bam METRICS_FILE=metrics_MQU_male_E8063.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8063.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8064\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037324_E8064_S102_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037324_E8064_S102_R2_filtered.fastq > MQU_male_E8064.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8064.sam MODE=SUMMARY O=MQU_male_E8064.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8064.sam OUTPUT=sorted_MQU_male_E8064.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8064.bam OUTPUT=dedup_MQU_male_E8064.bam METRICS_FILE=metrics_MQU_male_E8064.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8064.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8065\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037325_E8065_S103_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037325_E8065_S103_R2_filtered.fastq > MQU_male_E8065.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8065.sam MODE=SUMMARY O=MQU_male_E8065.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8065.sam OUTPUT=sorted_MQU_male_E8065.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8065.bam OUTPUT=dedup_MQU_male_E8065.bam METRICS_FILE=metrics_MQU_male_E8065.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8065.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8066\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037326_E8066_S104_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037326_E8066_S104_R2_filtered.fastq > MQU_male_E8066.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8066.sam MODE=SUMMARY O=MQU_male_E8066.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8066.sam OUTPUT=sorted_MQU_male_E8066.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8066.bam OUTPUT=dedup_MQU_male_E8066.bam METRICS_FILE=metrics_MQU_male_E8066.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8066.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8429\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037327_E8429_S91_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037327_E8429_S91_R2_filtered.fastq > MQU_male_E8429.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8429.sam MODE=SUMMARY O=MQU_male_E8429.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8429.sam OUTPUT=sorted_MQU_male_E8429.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8429.bam OUTPUT=dedup_MQU_male_E8429.bam METRICS_FILE=metrics_MQU_male_E8429.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8429.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8430\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037328_E8430_S100_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037328_E8430_S100_R2_filtered.fastq > MQU_male_E8430.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8430.sam MODE=SUMMARY O=MQU_male_E8430.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8430.sam OUTPUT=sorted_MQU_male_E8430.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8430.bam OUTPUT=dedup_MQU_male_E8430.bam METRICS_FILE=metrics_MQU_male_E8430.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8430.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8431\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037329_E8431_S101_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037329_E8431_S101_R2_filtered.fastq > MQU_male_E8431.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8431.sam MODE=SUMMARY O=MQU_male_E8431.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8431.sam OUTPUT=sorted_MQU_male_E8431.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8431.bam OUTPUT=dedup_MQU_male_E8431.bam METRICS_FILE=metrics_MQU_male_E8431.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8431.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8434\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037330_E8434_S94_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037330_E8434_S94_R2_filtered.fastq > MQU_male_E8434.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8434.sam MODE=SUMMARY O=MQU_male_E8434.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8434.sam OUTPUT=sorted_MQU_male_E8434.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8434.bam OUTPUT=dedup_MQU_male_E8434.bam METRICS_FILE=metrics_MQU_male_E8434.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8434.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8436\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037331_E8436_S102_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037331_E8436_S102_R2_filtered.fastq > MQU_male_E8436.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8436.sam MODE=SUMMARY O=MQU_male_E8436.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8436.sam OUTPUT=sorted_MQU_male_E8436.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8436.bam OUTPUT=dedup_MQU_male_E8436.bam METRICS_FILE=metrics_MQU_male_E8436.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8436.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8438\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037332_E8438_S96_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037332_E8438_S96_R2_filtered.fastq > MQU_male_E8438.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8438.sam MODE=SUMMARY O=MQU_male_E8438.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8438.sam OUTPUT=sorted_MQU_male_E8438.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8438.bam OUTPUT=dedup_MQU_male_E8438.bam METRICS_FILE=metrics_MQU_male_E8438.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8438.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8439\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037333_E8439_S103_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037333_E8439_S103_R2_filtered.fastq > MQU_male_E8439.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8439.sam MODE=SUMMARY O=MQU_male_E8439.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8439.sam OUTPUT=sorted_MQU_male_E8439.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8439.bam OUTPUT=dedup_MQU_male_E8439.bam METRICS_FILE=metrics_MQU_male_E8439.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8439.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8440\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037334_E8440_S104_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037334_E8440_S104_R2_filtered.fastq > MQU_male_E8440.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8440.sam MODE=SUMMARY O=MQU_male_E8440.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8440.sam OUTPUT=sorted_MQU_male_E8440.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8440.bam OUTPUT=dedup_MQU_male_E8440.bam METRICS_FILE=metrics_MQU_male_E8440.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8440.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8441\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037335_E8441_S105_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037335_E8441_S105_R2_filtered.fastq > MQU_male_E8441.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8441.sam MODE=SUMMARY O=MQU_male_E8441.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8441.sam OUTPUT=sorted_MQU_male_E8441.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8441.bam OUTPUT=dedup_MQU_male_E8441.bam METRICS_FILE=metrics_MQU_male_E8441.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8441.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8442\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037336_E8442_S100_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037336_E8442_S100_R2_filtered.fastq > MQU_male_E8442.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8442.sam MODE=SUMMARY O=MQU_male_E8442.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8442.sam OUTPUT=sorted_MQU_male_E8442.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8442.bam OUTPUT=dedup_MQU_male_E8442.bam METRICS_FILE=metrics_MQU_male_E8442.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8442.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8443\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037337_E8443_S7_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037337_E8443_S7_R2_filtered.fastq > MQU_male_E8443.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8443.sam MODE=SUMMARY O=MQU_male_E8443.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8443.sam OUTPUT=sorted_MQU_male_E8443.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8443.bam OUTPUT=dedup_MQU_male_E8443.bam METRICS_FILE=metrics_MQU_male_E8443.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8443.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8444\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037338_E8444_S102_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037338_E8444_S102_R2_filtered.fastq > MQU_male_E8444.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8444.sam MODE=SUMMARY O=MQU_male_E8444.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8444.sam OUTPUT=sorted_MQU_male_E8444.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8444.bam OUTPUT=dedup_MQU_male_E8444.bam METRICS_FILE=metrics_MQU_male_E8444.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8444.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8445\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037339_E8445_S103_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037339_E8445_S103_R2_filtered.fastq > MQU_male_E8445.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8445.sam MODE=SUMMARY O=MQU_male_E8445.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8445.sam OUTPUT=sorted_MQU_male_E8445.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8445.bam OUTPUT=dedup_MQU_male_E8445.bam METRICS_FILE=metrics_MQU_male_E8445.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8445.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8447\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037340_E8447_S104_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037340_E8447_S104_R2_filtered.fastq > MQU_male_E8447.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8447.sam MODE=SUMMARY O=MQU_male_E8447.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8447.sam OUTPUT=sorted_MQU_male_E8447.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8447.bam OUTPUT=dedup_MQU_male_E8447.bam METRICS_FILE=metrics_MQU_male_E8447.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8447.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8448\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037341_E8448_S105_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037341_E8448_S105_R2_filtered.fastq > MQU_male_E8448.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8448.sam MODE=SUMMARY O=MQU_male_E8448.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8448.sam OUTPUT=sorted_MQU_male_E8448.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8448.bam OUTPUT=dedup_MQU_male_E8448.bam METRICS_FILE=metrics_MQU_male_E8448.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8448.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8450\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037342_E8450_S106_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037342_E8450_S106_R2_filtered.fastq > MQU_male_E8450.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8450.sam MODE=SUMMARY O=MQU_male_E8450.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8450.sam OUTPUT=sorted_MQU_male_E8450.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8450.bam OUTPUT=dedup_MQU_male_E8450.bam METRICS_FILE=metrics_MQU_male_E8450.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8450.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8451\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037343_E8451_S107_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037343_E8451_S107_R2_filtered.fastq > MQU_male_E8451.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8451.sam MODE=SUMMARY O=MQU_male_E8451.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8451.sam OUTPUT=sorted_MQU_male_E8451.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8451.bam OUTPUT=dedup_MQU_male_E8451.bam METRICS_FILE=metrics_MQU_male_E8451.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8451.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8452\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037344_E8452_S108_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037344_E8452_S108_R2_filtered.fastq > MQU_male_E8452.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8452.sam MODE=SUMMARY O=MQU_male_E8452.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8452.sam OUTPUT=sorted_MQU_male_E8452.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8452.bam OUTPUT=dedup_MQU_male_E8452.bam METRICS_FILE=metrics_MQU_male_E8452.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8452.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8453\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037345_E8453_S109_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037345_E8453_S109_R2_filtered.fastq > MQU_male_E8453.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8453.sam MODE=SUMMARY O=MQU_male_E8453.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8453.sam OUTPUT=sorted_MQU_male_E8453.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8453.bam OUTPUT=dedup_MQU_male_E8453.bam METRICS_FILE=metrics_MQU_male_E8453.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8453.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8454\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037346_E8454_S110_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037346_E8454_S110_R2_filtered.fastq > MQU_male_E8454.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8454.sam MODE=SUMMARY O=MQU_male_E8454.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8454.sam OUTPUT=sorted_MQU_male_E8454.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8454.bam OUTPUT=dedup_MQU_male_E8454.bam METRICS_FILE=metrics_MQU_male_E8454.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8454.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8455\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037347_E8455_S111_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037347_E8455_S111_R2_filtered.fastq > MQU_male_E8455.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8455.sam MODE=SUMMARY O=MQU_male_E8455.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8455.sam OUTPUT=sorted_MQU_male_E8455.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8455.bam OUTPUT=dedup_MQU_male_E8455.bam METRICS_FILE=metrics_MQU_male_E8455.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8455.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8740\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037348_E8740_S112_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037348_E8740_S112_R2_filtered.fastq > MQU_male_E8740.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8740.sam MODE=SUMMARY O=MQU_male_E8740.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8740.sam OUTPUT=sorted_MQU_male_E8740.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8740.bam OUTPUT=dedup_MQU_male_E8740.bam METRICS_FILE=metrics_MQU_male_E8740.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8740.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8741\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037349_E8741_S113_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037349_E8741_S113_R2_filtered.fastq > MQU_male_E8741.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8741.sam MODE=SUMMARY O=MQU_male_E8741.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8741.sam OUTPUT=sorted_MQU_male_E8741.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8741.bam OUTPUT=dedup_MQU_male_E8741.bam METRICS_FILE=metrics_MQU_male_E8741.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8741.bam

bwa mem -t 40 -M -R "@RG\tID:group1\tSM:E8748\tPL:illumina\tLB:lib1\tPU:unit1" ../MQU_male.min500.fa \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037356_E8748_S120_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/2019/reads/run2/037356_E8748_S120_R2_filtered.fastq > MQU_male_E8748.sam
PicardCommandLine ValidateSamFile I=MQU_male_E8748.sam MODE=SUMMARY O=MQU_male_E8748.sam.txt
PicardCommandLine SortSam INPUT=MQU_male_E8748.sam OUTPUT=sorted_MQU_male_E8748.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_MQU_male_E8748.bam OUTPUT=dedup_MQU_male_E8748.bam METRICS_FILE=metrics_MQU_male_E8748.txt
PicardCommandLine BuildBamIndex INPUT=dedup_MQU_male_E8748.bam