!#/bin/bash -l

# This is the sicrp to Prepare a reference for use with BWA and GATK
# first we have to install all the software in the cluster

#BWA

cd ~

git clone https://github.com/lh3/bwa.git

cd bwa

make


#SAMTOOLS

git clone https://github.com/samtools/samtools.git

autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make
make install

#Picard


#Download of the complete genome of A. mexicanus to the cluster. ####TO:DO####




# Then we need to declare some enviromental variables in order to get simple paths

GENO=/Path/to/genome_reference.fa
BWA=/Path/to/BWA.version
SAM=/Path/to/SAMTOOLS
PIC=/Path/to/Picard.java

###############
###############
###############

#Generating the BWA index with BWA aligner

BWA index -a bwtsw GENO  																		#where -a bwtsw specifies that 
																								#we want to use the indexing algorithm 
																								#that is capable of handling very big
																								#genomes.

#Generation of the index fasta file with SAMTOOLS

SAM faidx GENO 																					#This creates a '.fa.fai' file 
																								#reference.fa.fai, with one record 
																								#per line for each of the contigs in 
																								#the FASTA reference file. Each record 
																								#is composed of the contig name, size, 
																								#location, basesPerLine and bytesPerLine.

#Generate the sequence dictionary

java -jar PIC CreateSequenceDictionary -R $GENO -O ../../Complete_genome/mexicanus.dict			#This creates a file called 
																								#reference.dict formatted like a SAM 
																								#header, describing the contents of 
																								#your reference FASTA file.
