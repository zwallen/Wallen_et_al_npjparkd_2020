#Get sequence count before, and after running DADA2

library(ShortRead); packageVersion("ShortRead")
library(dada2); packageVersion("dada2")

seqtab <- readRDS("../../Script_Output/Dataset2_Output/seqtab.rds")

directory.path <- "../../Sequences/Dataset2"

directories <- list.files(directory.path, pattern="5176-HP-Pool", full.names=TRUE)

seq.counts <- data.frame()

for (i in 1:length(directories)) {

	fastq.path <- directories[[i]]

	fastqFs <- list.files(fastq.path, pattern="_R1_001.fastq.gz", full.names=TRUE)
	fastqRs <- list.files(fastq.path, pattern="_R2_001.fastq.gz", full.names=TRUE)

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
}

write.table(seq.counts, "../../Script_Output/Dataset2_Output/SeqRemainingFrac.txt", sep="\t", row.names=F, quote=F)

