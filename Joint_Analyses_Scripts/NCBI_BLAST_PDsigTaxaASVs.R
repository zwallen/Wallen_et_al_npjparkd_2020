#BLAST ASVs for each PD-significant genus against NCBI's 16SMicrobial database
#Results reported in manuscript text and used in PubMed searches
date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")
library(seqinr); packageVersion("seqinr")
library(openxlsx); packageVersion("openxlsx")

#Read in dataset1 phyloseq object, remove phylogenetic tree
dataset1.ps <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")
dataset1.ps <- phyloseq(otu_table(dataset1.ps), tax_table(dataset1.ps), sample_data(dataset1.ps))
cat("\n","dataset1 phyloseq object summary:","\n")
dataset1.ps

#Read in dataset2 phyloseq object, remove phylogenetic tree
dataset2.ps <- readRDS("../PhyloseqObjects/Dataset2/phyloseq.rds")
dataset2.ps <- phyloseq(otu_table(dataset2.ps), tax_table(dataset2.ps), sample_data(dataset2.ps))
cat("\n","dataset2 phyloseq object summary:","\n")
dataset2.ps

#Create combined phyloseq object
ps <- merge_phyloseq(dataset1.ps, dataset2.ps)
cat("\n","Joint phyloseq object summary:","\n")
ps

#Transform to relative abundances 
ps.ra <- transform_sample_counts(ps, function(x){x/sum(x)})

#Define target taxa
target.taxa <- c("Porphyromonas", "Prevotella", "Corynebacterium_1", "Lactobacillus", "Bifidobacterium", "Agathobacter", "Lachnospiraceae_UCG-004", "Lachnospiraceae_ND3007_group", "Faecalibacterium", "Roseburia", "Blautia", "Oscillospira", "Fusicatenibacter", "Butyricicoccus", "Lachnospira")

#Run for ASVs within each target genus
wb <- createWorkbook()
for (i in 1:length(target.taxa)){

	#Subset phyloseq object for target Genus
	ps.ra.filt <- subset_taxa(ps.ra, Genus == target.taxa[i])
	cat("\n", "Running for ASVs belonging to", target.taxa[i], "\n")
	print(ps.ra.filt)

	#Get average proportions of each ASV
	results <- data.frame()
	for (ii in 1:length(taxa_names(ps.ra.filt))){
 
		avg.case.prop <- mean(otu_table(subset_samples(ps.ra.filt, case_control == "Case"))[,ii])
  		avg.control.prop <- mean(otu_table(subset_samples(ps.ra.filt, case_control == "Control"))[,ii])
  		ASV <- taxa_names(ps.ra.filt)[ii]
		
 		results <- rbind(results, data.frame(ASV = ASV, 
                                       Case_MRA = avg.case.prop, 
                                       Control_MRA = avg.control.prop))
	}
	#Calculate percent of genus an ASV makes up
	results$Perc_Genus_Case <- results$Case_MRA/sum(results$Case_MRA)
	results$Perc_Genus_Cont <- results$Control_MRA/sum(results$Control_MRA)
	
	#Add taxa designations to table
	full.results <- dplyr::inner_join(results, rownames_to_column(data.frame(tax_table(ps.ra.filt)), "ASV"), by = "ASV")

	#BLAST ASV sequences against NCBI 16S and pick out the top hit(s) for E-value, then % identity
	asv.vec <- c()
	for (j in 1:dim(full.results)[1]){
		
		#Create fasta file, blast against db, read back in, create column with formatted designations
		write.fasta(as.list(full.results$ASV)[[j]], names=paste("ASV",j,sep=""), file.out="ASV.fa")
		system("blastn -query ASV.fa -db ../Support_Files/ncbi-blast-2.9.0+/16SMicrobial -task megablast -outfmt \"6 stitle pident evalue\" > ASV_results.txt")
		asv.results <- read.table("ASV_results.txt", sep='\t', comment.char="", quote="\"")
		asv.results <- cbind(asv.results, paste(sapply(asv.results[,1], function(x){paste(strsplit(x, " ")[[1]][1], strsplit(x, " ")[[1]][2], sep=" ")}), " ", "(", asv.results[,2], ")", sep=""))
		
		#Filter for top hits based on e-value
		asv.top <- asv.results[asv.results[,3] == min(asv.results[,3]),]
		
		#Filter for top hits based on percent identity, remove duplicate entries if they exist, or concatenate different ones
		if (length(unique(asv.top[asv.top[,2] == max(asv.top[,2]), 4])) == 1){
			asv.vec <- append(asv.vec, unique(asv.top[asv.top[,2] == max(asv.top[,2]), 4]))
		}else{
			asv.vec <- append(asv.vec, paste(unique(asv.top[asv.top[,2] == max(asv.top[,2]), 4]), collapse="; "))
		}
	}

	#Add BLAST results to full results
	full.results$NCBI_Species <- asv.vec

	#Create worksheet in excel for results
	addWorksheet(wb, target.taxa[i])
	writeData(wb, target.taxa[i], full.results, keepNA=T)
}

#Write out results
saveWorkbook(wb, "../Script_Output/Joint_Analyses_Output/NCBI_BLAST_PDsigTaxaASVs.xlsx", overwrite=T)

#Clean up
system("rm ASV.fa")
system("rm ASV_results.txt")

