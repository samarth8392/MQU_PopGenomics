
##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     het.sh                       					###
###   					Estimate heterozygosity 						### 
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load samtools
module load R/3.4.2

cat /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/lists/az_sub_genic.list | while read -r LINE
do
	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 100 \
	-GL 1 -doSaf 1 -fold 1 -minQ 20 -minMapQ 20 \
	-anc /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/MQU_male.min500.fa \
	-i /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/align_nomerge/bam/${LINE}_genic.nomerge.bam \
	-out /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/az/${LINE}_genic_nomerge

	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 100 \
	/scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/az/${LINE}_genic_nomerge.saf.idx \
	> /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/az/${LINE}.genic.nomerge.ml
done

cat /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/lists/tx_sub_genic.list | while read -r LINE
do
	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 100 \
	-GL 1 -doSaf 1 -fold 1 -minQ 20 -minMapQ 20 \
	-anc /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/MQU_male.min500.fa \
	-i /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/align_nomerge/bam/${LINE}_genic.nomerge.bam \
	-out /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/tx/${LINE}_genic_nomerge

	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 100 \
	/scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/tx/${LINE}_genic_nomerge.saf.idx \
	> /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/tx/${LINE}.genic.nomerge.ml
done

cat /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/lists/nm_sub_genic.list | while read -r LINE
do
	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 100 \
	-GL 1 -doSaf 1 -fold 1 -minQ 20 -minMapQ 20 \
	-anc /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/MQU_male.min500.fa \
	-i /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/align_nomerge/bam/${LINE}_genic.nomerge.bam \
	-out /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/nm/${LINE}_genic_nomerge

	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 100 \
	/scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/nm/${LINE}_genic_nomerge.saf.idx \
	> /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/angsd/nomerge/het/nm/${LINE}.genic.nomerge.ml
done


# Whole genomic #

cat /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/lists/nm_sub_genic.bamlist | while read -r LINE
do
	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
	-GL 1 -doSaf 1 -fold 1 -minQ 20 -minMapQ 20 \
	-anc /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/MQU_male.min500.fa \
	-i /scratch/snyder/m/mathur20/MQU/ch2_redo/align_nomerge/final/${LINE}.mqu.recal.final.bam \
	-out /scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/az/${LINE}

	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 100 \
	/scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/az/${LINE}.saf.idx \
	> /scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/az/${LINE}.genomic.ml
done

cat /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/lists/tx_sub_genic.list | while read -r LINE
do
	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
	-GL 1 -doSaf 1 -fold 1 -minQ 20 -minMapQ 20 \
	-anc /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/MQU_male.min500.fa \
	-i /scratch/snyder/m/mathur20/MQU/ch2_redo/align_nomerge/final/${LINE}.mqu.recal.final.bam \
	-out /scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/tx/${LINE}

	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 \
	/scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/tx/${LINE}.saf.idx \
	> /scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/tx/${LINE}.genomic.ml
done

cat /scratch/snyder/m/mathur20/MQU/ch2_redo/genic/lists/nm_sub_genic.list | while read -r LINE
do
	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 200 \
	-GL 1 -doSaf 1 -fold 1 -minQ 20 -minMapQ 20 \
	-anc /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/MQU_male.min500.fa \
	-i /scratch/snyder/m/mathur20/MQU/ch2_redo/align_nomerge/final/${LINE}.mqu.recal.final.bam \
	-out /scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/nm/${LINE}

	/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 200 \
	/scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/nm/${LINE}.saf.idx \
	> /scratch/snyder/m/mathur20/MQU/ch2_redo/angsd/het/nm/${LINE}.genomic.ml
done