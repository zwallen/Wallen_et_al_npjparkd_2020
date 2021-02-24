#Perform correlation between PD associated genera and levodopa dose
#Results reported in manuscript text
date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")

taxa <- c("Porphyromonas","Prevotella","Corynebacterium_1","Agathobacter","Lachnospira","Lachnospiraceae_UCG-004","Lachnospiraceae_ND3007_group","Faecalibacterium","Butyricicoccus","Blautia","Roseburia","Fusicatenibacter","Oscillospira","Lactobacillus","Bifidobacterium")

phyloseq.object <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(phyloseq.object, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Transform counts to relative abundance
ps.genus.ra <- transform_sample_counts(ps.genus, function(x){x/sum(x)})
cat("\n","Genera counts transformed to proportions!","\n")

#extract target taxa
ps.genus.filtered <- subset_taxa(ps.genus.ra, Genus %in% taxa)
cat("\n","Extracting target taxa:","\n")
paste(tax_table(ps.genus.filtered)[,6])

#Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps.genus.filtered))
taxa.df <- rownames_to_column(taxa.df, "SampleID")

#Extract sample metadata from phyloseq object
meta.df <- data.frame(sample_data(ps.genus.filtered))

#Merge taxa.df and meta.df
taxa.meta.df <- dplyr::inner_join(taxa.df, meta.df, by = "SampleID")

### Loop Spearman's correlation for all taxa ###
results <- data.frame()

for (i in 2:length(taxa.df[1,])){
  
  OTU <- colnames(taxa.meta.df)[i]
 
  corr <- cor.test(taxa.meta.df[,i], taxa.meta.df$total_levodopa_dose, method = "spearman")
  
  results <- rbind(results, data.frame(Representative_ASV = OTU, 
                                       rho = corr$estimate,
                                       P = corr$p.value))
}

#Add taxa designations to table
results <- dplyr::inner_join(results, rownames_to_column(data.frame(tax_table(ps.genus.filtered)), "Representative_ASV"), by = "Representative_ASV")
write.table(results, "../Script_Output/Dataset1_Output/Corr_Ldopa_PDsigTaxa.txt", row.names = F, quote = F, sep = "\t")

