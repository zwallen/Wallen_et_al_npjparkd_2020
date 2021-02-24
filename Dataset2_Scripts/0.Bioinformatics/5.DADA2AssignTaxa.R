#Assign taxonomy to ASVs
library(dada2); packageVersion("dada2")
library(dplyr); packageVersion("dplyr")

#Load ASV table
seqtab <- readRDS("../../Script_Output/Dataset2_Output/seqtab.rds")

cat('\n',"Sequence variant table uploaded.", "\n")

#assign genus level taxonomy
taxa <- assignTaxonomy(seqtab, "../../Support_Files/silva_nr_v132_train_set.fa.gz", minBoot=80, multithread=TRUE)
cat("\n","Genus level taxonomy assigned.", "\n")

#add species level taxonomy
taxa <- addSpecies(taxa, "../../Support_Files/silva_species_assignment_v132.fa.gz", allowMultiple = TRUE)
cat("\n","Species level taxonomy added.", "\n")

#write taxa table to disk
saveRDS(taxa, "../../Script_Output/Dataset2_Output/ASVTaxaSILVA132minBoot80.rds")

#write taxa table
taxa.df <- add_rownames(as.data.frame(taxa), "ASV")
write.table(taxa.df, "../../Script_Output/Dataset2_Output/ASVTaxaSILVA132minBoot80.txt", sep="\t", row.names=F, quote=F)

cat("\n","Finished.", "\n")

