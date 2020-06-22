##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     smc.sh                        					###
###   			Demographic history estimation							### 
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N smc


while read -a line 
do
	for i in E6537 E6599 E6759 E6803 E6963 E7177 E7757
	do
		smc++ vcf2smc --cores 24 --ignore-missing --drop-first-last \
		-d $i $i -c 1000000 \
		/scratch/snyder/m/mathur20/MQU/2020/smc/vcf/az_chick_s2.vcf.gz \
		/scratch/snyder/m/mathur20/MQU/2020/smc/out_chick/AZ/out_chick.$i.${line[0]}.txt \
		${line[0]} \
		AZ:E6537,E6599,E6759,E6803,E6963,E7177,E7757
	done
done < /scratch/snyder/m/mathur20/MQU/2020/smc/chick_chr.txt


while read -a line 
do
	for i in E8430 E8431 E8438 E8444 E8445 E8447 E8740
	do
		smc++ vcf2smc --cores 24 --ignore-missing --drop-first-last \
		-d $i $i -c 1000000 \
		/scratch/snyder/m/mathur20/MQU/2020/smc/vcf/tx_chick_s2.vcf.gz \
		/scratch/snyder/m/mathur20/MQU/2020/smc/out_chick/TX/out_chick.$i.${line[0]}.txt \
		${line[0]} \
		TX:E8430,E8431,E8438,E8444,E8445,E8447,E8740
	done
done < /scratch/snyder/m/mathur20/MQU/2020/smc/chick_chr.txt

while read -a line 
do
	for i in E8451 E8452 E8453 E8454 E8455 E8741 E8748
	do
		smc++ vcf2smc --cores 24 --ignore-missing --drop-first-last \
		-d $i $i -c 1000000 \
		/scratch/snyder/m/mathur20/MQU/2020/smc/vcf/nm_chick_s2.vcf.gz \
		/scratch/snyder/m/mathur20/MQU/2020/smc/out_chick/NM/out_chick.$i.${line[0]}.txt \
		${line[0]} \
		NM:E8451,E8452,E8453,E8454,E8455,E8741,E8748
	done
done < /scratch/snyder/m/mathur20/MQU/2020/smc/chick_chr.txt


#for p in AZ TX NM
do
	smc++ cv --cores 400 --base ${p}_th_1300_rp_6 -v \
	-o /scratch/snyder/m/mathur20/MQU/2020/smc/estimate_out/chick/r6/ \
	--thinning 1300 --em-iterations 1000 \
	--timepoints 1e3 1e6 \
	--algorithm Powell --multi --ftol 1e-7 --xtol 1e-7 \
	--regularization-penalty 6 \
	--unfold --knots 10 \
	--spline piecewise --folds 5 \
	3.13597726022324e-09  \
	/scratch/snyder/m/mathur20/MQU/2020/smc/out_chick/${p}/*.txt
done

smc++ plot /scratch/snyder/m/mathur20/MQU/2020/smc/estimate_out/chick/r6/plots/AZNMTX_th_1300_rp_6.pdf \
-g 1.5 --linear -c \
/scratch/snyder/m/mathur20/MQU/2020/smc/estimate_out/chick/r6/fold0/AZ_th_1300_rp_6.final.json \
/scratch/snyder/m/mathur20/MQU/2020/smc/estimate_out/chick/r6/fold0/TX_th_1300_rp_6.final.json \
/scratch/snyder/m/mathur20/MQU/2020/smc/estimate_out/chick/r6/fold0/NM_th_1300_rp_6.final.json 