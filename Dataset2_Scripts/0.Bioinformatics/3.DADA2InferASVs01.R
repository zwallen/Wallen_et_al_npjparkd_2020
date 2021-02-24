#Infer Sequence Variants by learning the error rates, depreplicating sequences to unique sequences, and performing dada2 algorithm

library(dada2); packageVersion("dada2")
getN <- function(x) sum(getUniques(x))

#File parsing

filtpath <- "../../Script_Output/Dataset2_Output/PrimerTrimmedSeqs/QualityFiltSeqs"

filtFs <- list.files(filtpath, pattern="_01_R1_001.fastq.gz", full.names=TRUE)
filtRs <- list.files(filtpath, pattern="_01_R2_001.fastq.gz", full.names=TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

if (!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match")

names(filtFs) <- sample.names
names(filtRs) <- sample.namesR

#Learn forward read error rates
set.seed(1234)
errF <- learnErrors(filtFs, nbases=Inf, multithread=TRUE)
#Learn reverse read error rates
set.seed(1234)
errR <- learnErrors(filtRs, nbases=Inf, multithread=TRUE)

#Make pdfs of error rate plots
#Forward
pdf("../../Script_Output/Dataset2_Output/ErrorRatePlotF_01.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()
#Reverse
pdf("../../Script_Output/Dataset2_Output/ErrorRatePlotR_01.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

#Sequence variant inference and merging of paired end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sample in sample.names) {
	cat("Processing:", sample, "\n")
		derepF <- derepFastq(filtFs[[sample]])
		ddF <- dada(derepF, err=errF, multithread=TRUE)
		cat('Forward reads remaining after denoising:',getN(ddF),'\n')
		derepR <- derepFastq(filtRs[[sample]])
		ddR <- dada(derepR, err=errR, multithread=TRUE)
		cat('Reverse reads remaining after denoising:',getN(ddR),'\n')
		merger <- mergePairs(ddF, derepF, ddR, derepR)
		cat('Reads remaining after merging:',getN(merger),'\n')
		mergers[[sample]] <- merger
}
rm(derepF); rm(derepR)

#Construct initial sequence variant table
seqtab <- makeSequenceTable(mergers)

#Get sequence length distribution:
cat("\n","ASV sequence length distribution:", "\n")
table(nchar(getSequences(seqtab)))

#Extract only certain length of sequence at appropriate length (like cutting a gel band)
seqtab2 <- seqtab[,nchar(colnames(seqtab))%in%seq(250,256)]

saveRDS(seqtab2, "../../Script_Output/Dataset2_Output/seqtab01.rds")
