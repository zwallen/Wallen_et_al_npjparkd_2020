#!/bin/bash

#Check if output directory for primer trimmed sequences already exists
#If yes, delete old directory and make new one
#If no, make new directory
if [ -d "../../Script_Output/Dataset1_Output/PrimerTrimmedSeqs" ]
then
	rm -r ../../Script_Output/Dataset1_Output/PrimerTrimmedSeqs
	mkdir ../../Script_Output/Dataset1_Output/PrimerTrimmedSeqs
else
	mkdir ../../Script_Output/Dataset1_Output/PrimerTrimmedSeqs
fi

#Run cutadapt to trim primer sequences
for f in ../../Sequences/Dataset1/

do
	cd $f
	for i in *R1_001.fastq.gz
	do
		FILE_NAME=$(echo $i | awk -F '_' '{print $1}')

		cutadapt -g GTGCCAGCMGCCGCGGTAA \
		-G GGACTACHVGGGTWTCTAAT \
		-o ../../Script_Output/Dataset1_Output/PrimerTrimmedSeqs/${FILE_NAME}_R1_001.fastq.gz \
		-p ../../Script_Output/Dataset1_Output/PrimerTrimmedSeqs/${FILE_NAME}_R2_001.fastq.gz \
		${FILE_NAME}_R1_001.fastq.gz \
		${FILE_NAME}_R2_001.fastq.gz \
		--minimum-length 147 --maximum-length 151 \
		> ../../Script_Output/Dataset1_Output/PrimerTrimmedSeqs/${FILE_NAME}.log
	done
done
