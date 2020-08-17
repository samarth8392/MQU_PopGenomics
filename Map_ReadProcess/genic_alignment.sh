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

# GFF TO BED using bedops#

awk '{gsub(/;$/,"");print}' genome.all.gff > revisedgenome.all.gff 
cat revisedgenome.all.gff | awk '$3 =="gene" {print $0}' > revisedgenome.all_justgenes.gff 
gff2bed < revisedgenome.all_justgenes.gff > revisedgenome.all_justgenes.bed 


# Extract genic BAM using samtools view #
while read -a line
do
	samtools view -L /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/align/revisedgenome.all_justgenes.bed \
	-o /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/align_nomerge/bam/${line[0]}_genic.nomerge.bam \
	/scratch/snyder/m/mathur20/MQU/ch2_redo/align_nomerge/final/${line[0]}.mqu.recal.final.bam
done < /scratch/snyder/m/mathur20/MQU/ch2_redo/lists/all.list


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


