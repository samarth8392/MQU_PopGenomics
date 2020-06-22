###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###########################################################################
###########################################################################
###                     mitogenome_alignment.sh                        	### 
###																		###
###  To assemble mitogenomes from WGS data                              ###
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N mitogenome_alignment

cd $PBS_O_WORKDIR
module purge

# STEP 1 #
# Denovo assmebly of mitogenome using MitoBim
# (Already done for MQU; Mathur et al. 2019)


## STEP 2: Find numts

# To annotate nuclear mitochondrial translocations (NUMTs) in the jaguar reference genome, we used a BLAST-based approach. 
# Using the complete jaguar mitochondrial genome (generated in this study) as a query, matches were sought using the following 
# search parameters: 
#(i) a hit with at least 16 bp; 
#(ii) E value threshold of1 × 10−10; 
#(iii) noDUST filter query; 
#(iv) cost of0 to open a gap and 2 to extend it; 
#(v) X dropoff value of 40, for preliminary gapped extensions; and 
#(vi) reward for a matchand penalt for a mismatch of 1

module purge
module load bioinfo
module load blast/2.5.0+

dustmasker -in /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female.min500_renamed.fasta \
-infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out MQU_female_renamed_refdb.asnb

makeblastdb -in /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female.min500_renamed.fasta \
-input_type fasta -parse_seqids -dbtype nucl -out MQU_female_renamed_refdb 

blastn -db /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female_renamed_refdb  \
-evalue 0.0001 -outfmt '10 qseqid sseqid stitle pident evalue qcovs' \
-gapopen 0 -gapextend 2 -word_size 16 -xdrop_gap 40 \
-reward 1 -penalty -1 \
-query /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female_mito_renamed.fasta \
-out /scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/MQU_potential_numts.txt

# STEP 3: Get numt fasta sequences from the nuclear genome

module load bioinfo
module load biopython/3.6.5

/scratch/snyder/m/mathur20/MQU/2019/jobcodes/faSomeRecords \
/scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female.min500_renamed.fasta \
/scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/MQU_genome_numt_id.txt \
/scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/MQU_numts.fasta

## STEP 4: Extract mito-specific reads that map to numts sequences

cd /scratch/snyder/m/mathur20/MQU/2019/mito/round2/mito_reads/step1/
for prefix in $(ls *.fq | sed -r 's/_step1_mitoMQU_female_mito_[12][.]fq//'  | uniq)
do
	cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_jobs/
	echo 	"#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N ${prefix}_Step4
cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load blast/2.5.0+
module load MITObim/1.8
module load BBMap/37.93

/group/bioinfo/apps/apps/MITObim-1.8/misc_scripts/interleave-fastqgz-MITOBIM.py \
/scratch/snyder/m/mathur20/MQU/2019/mito/round2/mito_reads/step1/${prefix}_step1_mitoMQU_female_mito_1.fq \
/scratch/snyder/m/mathur20/MQU/2019/mito/round2/mito_reads/step1/${prefix}_step1_mitoMQU_female_mito_2.fq > /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/${prefix}_interleaved_wnumts.fq


bbsplit.sh -Xmx4g build=1 \
ref_MQU=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/MQU_numts.fasta \
in=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/${prefix}_interleaved_wnumts.fq \
out_MQU=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/${prefix}_only_numts.fastq" > ${prefix}_step4.job
done

## STEP 5: Remove sequences that match numts from mito-specific reads

cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/
for prefix in $(ls *.fq | sed -r 's/_interleaved_wnumts[.]fq//'  | uniq)
do
	cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_jobs/
	echo 	"#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N Step5_${prefix}
cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load blast/2.5.0+
module load BBMap/37.93

grep \"@\" /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/${prefix}_only_numts.fastq > /scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/${prefix}_numt_read_ids
cat /scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/${prefix}_numt_read_ids | tr -d '@' > /scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/${prefix}_numt_ids

filterbyname.sh in=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/${prefix}_interleaved_wnumts.fq \
out=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/${prefix}_interleaved_without_numts.fq \
names=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/${prefix}_numt_read_ids \
include=f ow=t substring=header ths=t" > ${prefix}_step5.job
done

## STEP 6: Update the mitogenome sequence
# Update the reference mitogenome  #

module purge
module load bioinfo
module load BBMap/37.93
module load MITObim/1.8

bbsplit.sh -Xmx4g build=1 \
ref_MQU=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/MQU_numts.fasta \
in=/scratch/snyder/m/mathur20/MQU/2018/genome_assembly/mitobim/female_mito/iteration113/MQU_female-readpool-it113.fastq \
out_MQU=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_mito_ref_only_numts.fastq
grep "@" /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_mito_ref_only_numts.fastq > /scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/ref_numt_read_ids
cat /scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/ref_numt_read_ids | tr -d '@' > /scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/ref_numt_ids

filterbyname.sh in=/scratch/snyder/m/mathur20/MQU/2018/genome_assembly/mitobim/female_mito/iteration113/MQU_female-readpool-it113.fastq \
out=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_mito_ref_without_numts.fq \
names=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/numts/ref_numt_ids \
include=f ow=t substring=header ths=t

mkdir /scratch/snyder/m/mathur20/MQU/2019/mito/round3/update_mito
cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/update_mito

MITObim.pl -start 1 -end 500 -ref MQU_update_mito \
-sample MQU_update_mito --clean \
--quick /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female_mito.fasta \
-readpool /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_mito_ref_without_numts.fastq &> log

## Step 7: Map all individuals to updated mitogenome

cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/no_numts/
for prefix in $(ls| sed -r 's/_mito_reads//'  | uniq)
do
	cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_jobs/
	echo 	"#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N Step7_${prefix}
cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load bwa
module load BBMap
module load picard-tools
module load samtools
module load BBTools/35.69

#rename.sh -Xmx4g ow=t in=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/no_numts/${prefix}_mito_reads/${prefix}_interleaved_without_numts.fq \
#out=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/no_numts/${prefix}_mito_reads/${prefix}_interleaved_without_numts_renamed.fq prefix=${prefix}

#bbnorm.sh -Xmx40g threads=4 in=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/no_numts/${prefix}_mito_reads/${prefix}_interleaved_without_numts_renamed.fq \
#out=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/no_numts/${prefix}_mito_reads/${prefix}_interleaved_without_numts_renamed_mindepth50.fq \
#target=99999999 mindepth=50 passes=2

#bwa index /scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female_mito_nonumts_trim.fasta

bwa mem -t 4 -M -R \"@RG\tID:group1\tSM:${prefix}\tPL:illumina\tLB:lib1\tPU:unit1\" \
/scratch/snyder/m/mathur20/MQU/2019/mito/round3/MQU_ref/MQU_female_mito_nonumts_trim.fasta \
/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_reads/no_numts/${prefix}_mito_reads/${prefix}_interleaved_without_numts_renamed.fq > /scratch/snyder/m/mathur20/MQU/2019/mito/round3/alignments/${prefix}_MQU_nonumt_mito.sam

cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/alignments/

PicardCommandLine ValidateSamFile I=${prefix}_MQU_nonumt_mito.sam MODE=SUMMARY O=${prefix}_MQU_nonumt_mito.sam_samfile.txt
PicardCommandLine SortSam INPUT=${prefix}_MQU_nonumt_mito.sam OUTPUT=sorted_${prefix}_MQU_nonumt_mito.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_${prefix}_MQU_nonumt_mito.bam OUTPUT=dedup_sorted_${prefix}_MQU_nonumt_mito.bam METRICS_FILE=metrics_${prefix}_mito
PicardCommandLine BuildBamIndex INPUT=dedup_sorted_${prefix}_MQU_nonumt_mito.bam
samtools depth -a dedup_sorted_${prefix}_MQU_nonumt_mito.bam > ${prefix}_mito_coverage.txt"> ${prefix}_step7.job
done


################## RUN JOB CODES ############################


cd /scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_jobs/
for prefix in $(ls *.job | sed -r 's/_step[012345678][.]job//' | uniq)
do
	qsub ${prefix}_step7.job
done
