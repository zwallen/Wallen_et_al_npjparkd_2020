#Detect Genera significantly different between PD patients and controls using Kruskal-Wallis test
#No adjustment for covariates; normalizing taxa counts by total library sizes
#Results reported in Table 2 and Supplementary Table 4
date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")

ps <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(ps, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Transform counts to relative abundance
ps.genus.ra <- transform_sample_counts(ps.genus, function(x){x/sum(x)})
cat("\n","Genera counts transformed to proportions!","\n")

#Filter out any unclassified groups
ps.genus.noNA <- subset_taxa(ps.genus.ra, !is.na(Genus))
cat("\n","Removing unclassified groups:", "\n")
ps.genus.noNA

#Filter taxa found in less than 10% of samples
ps.genus.filtered <- filter_taxa(ps.genus.noNA, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)
cat("\n","Extracting Genera found in at least 10% of samples:","\n")
ps.genus.filtered

#Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps.genus.filtered))
taxa.df <- rownames_to_column(taxa.df, "SampleID")

#Extract sample metadata from phyloseq object
meta.df <- data.frame(sample_data(ps.genus.filtered))
cat("\n","Case/control included in analysis:","\n")
table(meta.df$case_control)

#Merge taxa.df and meta.df
taxa.meta.df <- dplyr::inner_join(taxa.df, meta.df, by = "SampleID")

### Loop Kruskal-Wallis test for all taxa ###
results <- data.frame()

for (i in 2:length(taxa.df[1,])){
 
  avg.case.prop <- mean(subset(taxa.meta.df, case_control == "Case")[,i])
  avg.control.prop <- mean(subset(taxa.meta.df, case_control == "Control")[,i])
  abundance.ratio <- avg.case.prop/avg.control.prop
  OTU <- colnames(taxa.meta.df)[i]
 
  kw.pvalue <- kruskal.test(taxa.meta.df[,i] ~ as.factor(case_control), 
                            data = taxa.meta.df)$p.value
  results <- rbind(results, data.frame(Representative_ASV = OTU, 
                                       PD_MRA = avg.case.prop, 
                                       Control_MRA = avg.control.prop, 
                                       FC = abundance.ratio, 
                                       P = kw.pvalue))
}
#Perform Benjamini-Hochberg FDR correction for KW pvalues
results$FDR_BH <- p.adjust(results$P, method = 'BH')

#Add taxa designations to table
results <- dplyr::inner_join(results, rownames_to_column(data.frame(tax_table(ps.genus.filtered)), "Representative_ASV"), by = "Representative_ASV")
write.table(results, "../Script_Output/Dataset1_Output/KW_MWAS.txt", quote = F, sep="\t", row.names = F)
