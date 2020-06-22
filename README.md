# MQU_PopGenomics
## Montezuma Quail Population Genomics

This repository contains scripts for the publication:

Samarth Mathur, J. Andrew DeWoody. Genetic load has potential in large populations but is realized in small populations. 
Authorea. May 14, 2020.DOI: 10.22541/au.158941448.85174067

## Folders
This repository contains the following scripts within each folder

### Map_ReadProcess
#### Scripts for aligning MQU reads to a reference genome and processing reads according to GATK Best Practices

##### alignment.sh 
Aligning reads to a reference assembly, sorting, marking dupicates, realinging around indels, and base quality score recalibration

##### genic_alignment.sh 
Extracting genic reads and mapping it to MQU genic regions (as annotated in Mathur et al. 2019)
