#!/bin/bash -l        
#PBS -l walltime=8:00:00,nodes=1:ppn=8,mem=10gb 
#PBS -m abe 
#PBS -M mgarduos@umn.edu 

#TRIMMING with TRIMMOMATIC.

#Download Trimmomatic from GitHug source

cd ~

mkdir trimmomatic

cd trimmomatic

git clone https://github.com/timflutre/trimmomatic.git

make

make check

make install

# stablish all my locations and make new folders

TRIM=PATH/TO/TRIMMOMATIC
NEXT=PATH/TO/NEXTERA_ADAPTERS/
DATA=PATH/TO/MY/DATA

mkdir ./TRIM

#Then i can make two things
#1) stablish a 'for' loop, in order to run all my data 

for i in $(ls *fastq.gz); do 
	java -jar $TRIM SE -phred33 $i ./TRIM/TRIM_$i ILLUMINACLIP:$NEXT:2:30:10 HEADCROP:14 SLIDINGWINDOW:4:15 MINLEN:75; 
done

#2) Stablis a parallel conection 
