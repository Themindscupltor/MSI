#!/bin/bash -l
#PBS -l walltime=8:00:00,nodes=1:ppn=8,mem=10gb 
#PBS -m abe 
#PBS -M mgarduos@umn.edu 

#With this script we can run On my personal computer

#Stablishing the genome reference 

GENO=/path/to/genome_reference.fastq							#This had to be indexed previously using BWA, SAMTOOLS and PICARD
	
#Using BWA MEM with my pilot samples.

BWA:path/to/BWA

#Move to my trimmed data

cd /to/my/trimmed/data
mkdir ALN

#Run BWA MEM

$BWA mem $GENO TRIM_Chiquitita2_S9_R1_001.fastq.gz > ./ALN/ALN_SE_CHI.sam

#Then i have to give some kind of format to my samples un order to get adapt to the SNP caller

PIC=path/to/picard.jar

java -jar $PIC SamFormatConverter I=ALN_SE_CHI.sam O=ALN_SE_CHI.bam

# Order the aligned reads by position using the reference ****SORTSAM****

java -jar $PIC SortSam I= ./ALN_SE_CHI.bam O= ./ALN_SE_CHI_sorted.bam SO=coordinate

#Add the ReadGroup Field to the .bam file using AddOrReplaceReadGroups command

java -jar picard.jar AddOrReplaceReadGroups I=ALN_SE_CHI_sorted.bam O=ALN_SE_CHI_sorted_RG.bam \
ID=muestra LB=Paired-end PL=Illumina PU=Unknown SM=muestra

#Indexing the aligment 

samtools index ALN_SE_CHI_sorted_RG.bam
