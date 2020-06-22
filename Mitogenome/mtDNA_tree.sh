###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 04/03/19                  Last Modified: 11/26/19 ###
###########################################################################
###########################################################################
###                     mtDNA_tree.sh                        			### 
### Create phylogenetic trees from mtDNA								###
###########################################################################

#!/bin/sh -l
#PBS -q gcore
#PBS -l nodes=1:ppn=4,naccesspolicy=singleuser
#PBS -l walltime=336:00:00
#PBS -N mtDNA_tree

cd $PBS_O_WORKDIR
module purge
module load bioinfo
module load clustalw/2.1

clustalw -ALIGN -TREE -BOOTSTRAP=1000 -TYPE=DNA -OUTPUT=NEXUS -OUTPUTTREE=nj \
-INFILE=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/fasta/all_MQU_filtered.fa \
-OUTFILE=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_tree/post_nov19/all_MQU_filtered

clustalw -ALIGN -TREE -BOOTSTRAP=1000 -TYPE=DNA -OUTPUT=NEXUS -OUTPUTTREE=nj \
-INFILE=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/fasta/only_AZ_filtered.fa \
-OUTFILE=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_tree/post_nov19/only_AZ_filtered

clustalw -ALIGN -TREE -BOOTSTRAP=1000 -TYPE=DNA -OUTPUT=NEXUS -OUTPUTTREE=nj \
-INFILE=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/fasta/only_TX_filtered.fa \
-OUTFILE=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_tree/post_nov19/only_TX_filtered

clustalw -ALIGN -TREE -BOOTSTRAP=1000 -TYPE=DNA -nexus_OUTPUT=NEXUS -OUTPUTTREE=nj \
-INFILE=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/fasta/only_NM_filtered.fa \
-OUTFILE=/scratch/snyder/m/mathur20/MQU/2019/mito/round3/mito_tree/post_nov19/only_NM_filtered