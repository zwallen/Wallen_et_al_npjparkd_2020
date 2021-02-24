#!/bin/bash

#Check if output directory for primer trimmed sequences already exists
#If yes, proceed with processing
#If no, give error
if [ -d "../../Script_Output/Dataset2_Output/PrimerTrimmedSeqs" ]
then
	:
else
	echo "Error: no directory exists to place primer trimmed sequences, did you run 1.RemovePrimer1.sh script first?"
	exit 1
fi

#Run cutadapt to trim primer sequences
for f in ../../Sequences/Dataset2/5176-HP-Pool_05*

do
	cd $f
	for i in *R1_001.fastq.gz
	do
		FILE_NAME=$(echo $i | awk -F '_' '{print $1}')

		cutadapt -g GTGCCAGCMGCCGCGGTAA \
		-G GGACTACHVGGGTWTCTAAT \
		-o ../../../Script_Output/Dataset2_Output/PrimerTrimmedSeqs/${FILE_NAME}_05_R1_001.fastq.gz \
		-p ../../../Script_Output/Dataset2_Output/PrimerTrimmedSeqs/${FILE_NAME}_05_R2_001.fastq.gz \
		${FILE_NAME}_R1_001.fastq.gz \
		${FILE_NAME}_R2_001.fastq.gz \
		--minimum-length 230 --maximum-length 233 \
		> ../../../Script_Output/Dataset2_Output/PrimerTrimmedSeqs/${FILE_NAME}.log
	done
done
