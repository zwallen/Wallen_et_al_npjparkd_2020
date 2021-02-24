#Get beginning sequence count, and sequence count after running DADA2

library(ShortRead); packageVersion("ShortRead")
library(dada2); packageVersion("dada2")

seqtab <- readRDS("../../Script_Output/Dataset1_Output/seqtab.rds")

directory.path <- "../../Sequences/Dataset1"

seq.counts <- data.frame()

#for (i in 1:length(directories)) {

	fastqFs <- list.files(directory.path, pattern="_R1_001.fastq.gz", full.names=TRUE)
	fastqRs <- list.files(directory.path, pattern="_R2_001.fastq.gz", full.names=TRUE)

	sample.names <- sapply(strsplit(basename(fastqFs), "_"), `[`, 1)
	sample.namesR <- sapply(strsplit(basename(fastqRs), "_"), `[`, 1)

	if (!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match")

	names(fastqFs) <- sample.names
	names(fastqRs) <- sample.namesR

	#Sequence tabulation 
	for(sample in sample.names) {
		cat("Processing:", sample, "\n")
			fastqF <- readFastq(fastqFs[[sample]])
			fastqR <- readFastq(fastqRs[[sample]])
			if (length(fastqF) != length(fastqR)) stop("Forward and reverse files do not have same sequence count")
			initial.seq.count <- length(fastqF)
			final.seq.count <- sum(seqtab[sample,])
			frac.seq.remain <- round((final.seq.count/initial.seq.count), 2)
			seq.counts <- rbind(seq.counts, data.frame(SampleID=sample, InitialCount=initial.seq.count, FinalCount=final.seq.count, FracSeqRemain=frac.seq.remain))
	}
#}

write.table(seq.counts, "../../Script_Output/Dataset1_Output/SeqRemainingFrac.txt", sep="\t", row.names=F, quote=F)

