#!/bin/bash
#Download sequences from SRA (BioProject PRJNA601994) using SRAToolkit

### Make directories to place sequences into ###

#Dataset 1
mkdir Dataset1

#Dataset 2
mkdir Dataset2
for i in {1..6}
do
  mkdir Dataset2/5176-HP-Pool_0${i}
done

### Download sequences from SRA into their respective directories ###

#Dataset 1
grep "dataset_1" SraRunTable.txt | \
while read line; do

    OUT_DIR=Dataset1
    SRR_ID=$(echo $line | grep -o '\S*SRR\S*' | uniq)
    SAMP_ID=$(echo $line | grep -o '\S*10122\S*' | uniq)

    prefetch $SRR_ID --output-directory $OUT_DIR
    fastq-dump ${OUT_DIR}/${SRR_ID}.sra --split-files --outdir $OUT_DIR --gzip

    mv ${OUT_DIR}/${SRR_ID}_1.fastq.gz ${OUT_DIR}/${SAMP_ID}_R1_001.fastq.gz
    mv ${OUT_DIR}/${SRR_ID}_2.fastq.gz ${OUT_DIR}/${SAMP_ID}_R2_001.fastq.gz

    rm ${OUT_DIR}/${SRR_ID}.sra
    rm ~/ncbi/public/sra/*.cache

done

#Dataset 2
for i in {1..6}
do

    grep "dataset_2_${i}" SraRunTable.txt | \
    while read line; do

        OUT_DIR=Dataset2/5176-HP-Pool_0${i}
        SRR_ID=$(echo $line | grep -o '\S*SRR\S*' | uniq)
        SAMP_ID=$(echo $line | grep -o '\S*5176\S*' | uniq)

        prefetch $SRR_ID --output-directory $OUT_DIR
        fastq-dump ${OUT_DIR}/${SRR_ID}.sra --split-files --outdir $OUT_DIR --gzip

        mv ${OUT_DIR}/${SRR_ID}_1.fastq.gz ${OUT_DIR}/${SAMP_ID}_R1_001.fastq.gz
        mv ${OUT_DIR}/${SRR_ID}_2.fastq.gz ${OUT_DIR}/${SAMP_ID}_R2_001.fastq.gz

        rm ${OUT_DIR}/${SRR_ID}.sra
        rm ~/ncbi/public/sra/*.cache

    done

done

