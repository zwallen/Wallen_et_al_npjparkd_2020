#Detect Genera significantly different between PD patients and controls using Kruskal-Wallis test
#This time with cluster 1 genera aggregated together (leaving PD associated genera out)
#Results reported in manuscript text
date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")

phyloseq.object <- readRDS("../PhyloseqObjects/Dataset2/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(phyloseq.object, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Extract genera in cluster 1
partial.clust <- c("Anaerococcus", "Campylobacter", "Ezakiella", "Fastidiosipila", "Finegoldia", "Lawsonella", "Mobiluncus", "Mogibacterium", "Murdochiella", "Negativicoccus", "Peptoniphilus", "Prevotella_6", "S5-A14a", "Varibaculum","Corynebacteriaceae_NA")
ps.genus.part.clust <- subset_taxa(ps.genus, Genus %in% partial.clust | (Family == "Corynebacteriaceae" & is.na(Genus)))

#Filter out cluster genera from main phyloseq object
ps.genus.sub.part.clust <- subset_taxa(ps.genus, !(colnames(otu_table(ps.genus)) %in% colnames(otu_table(ps.genus.part.clust))))

#Add cluster abundance as one group
ps.genus.part.clust <- tax_glom(ps.genus.part.clust, taxrank = "Kingdom", NArm = F)
taxa_names(ps.genus.part.clust) <- "Partial_Cluster1"
ps.genus.sub.part.clust <- phyloseq(otu_table(ps.genus.sub.part.clust), tax_table(ps.genus.sub.part.clust), sample_data(ps.genus.sub.part.clust))
ps.genus.sub.part.clust <- merge_phyloseq(ps.genus.sub.part.clust, ps.genus.part.clust)
cat("\n","Phyloseq object with partial cluster 1 added:","\n")
ps.genus.sub.part.clust

#Transform counts to relative abundance
ps.genus.sub.part.clust.ra <- transform_sample_counts(ps.genus.sub.part.clust, function(x){x/sum(x)})
cat("\n","Genera counts transformed to proportions!","\n")

#Filter out any unclassified groups
ps.genus.sub.part.clust.noNA <- subset_taxa(ps.genus.sub.part.clust.ra, !is.na(Genus) | taxa_names(ps.genus.sub.part.clust.ra) == "Partial_Cluster1")
cat("\n","Removing unclassified groups:", "\n")
ps.genus.sub.part.clust.noNA

#Filter taxa found in less than 10% of samples
ps.genus.sub.part.clust.filtered <- filter_taxa(ps.genus.sub.part.clust.noNA, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
cat("\n","Extracting Genera found in greater than 10% of samples:","\n")
ps.genus.sub.part.clust.filtered

#Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps.genus.sub.part.clust.filtered))
taxa.df <- rownames_to_column(taxa.df, "SampleID")

#Extract sample metadata from phyloseq object
meta.df <- data.frame(sample_data(ps.genus))
meta.df <- meta.df[,-1]
colnames(meta.df)[1] <- "SampleID"
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
results <- dplyr::inner_join(results, rownames_to_column(data.frame(tax_table(ps.genus.sub.part.clust.filtered)), "Representative_ASV"), by = "Representative_ASV")

write.table(results, "../Script_Output/Dataset2_Output/KW_MWAS_Cluster1.txt", quote = F, sep="\t", row.names = F)
