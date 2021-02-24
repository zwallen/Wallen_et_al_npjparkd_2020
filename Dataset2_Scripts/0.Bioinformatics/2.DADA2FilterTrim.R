#Perform filtering of paired-end forward and reverse reads using DADA2

library(dada2); packageVersion("dada2")

#Parse forward and reverse read files

path <- "../../Script_Output/Dataset2_Output/PrimerTrimmedSeqs"

filtpath <- file.path(path, "QualityFiltSeqs")

fastqFs <- sort(list.files(path, pattern="_R1_001.fastq.gz"))
fastqRs <- sort(list.files(path, pattern="_R2_001.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#Filtering of sequence files

filterAndTrim(fwd=file.path(path, fastqFs), filt=file.path(filtpath, fastqFs), rev=file.path(path, fastqRs), filt.rev=file.path(filtpath, fastqRs), truncLen=c(228,203), maxEE=2, truncQ=2, maxN=0, rm.phix=T, compress=T, verbose=T, multithread=F, n=1e5)


