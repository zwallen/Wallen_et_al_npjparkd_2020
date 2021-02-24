#Detect Genera significantly different between PD patients and controls using ANCOM adjusting for covariates
#This time with cluster 1 genera aggregated together (leaving PD associated genera out)
#Results reported in manuscript text
date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")
library(dplyr); packageVersion("dplyr")
source("../Support_Files/ANCOM_updated_code.R")

ps <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(ps, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Extract genera in cluster 1
partial.clust <- c("Anaerococcus", "Campylobacter", "Ezakiella", "Finegoldia", "Murdochiella", "Peptoniphilus", "Varibaculum")
ps.genus.part.clust <- subset_taxa(ps.genus, Genus %in% partial.clust)

#Filter out cluster genera from main phyloseq object
ps.genus.sub.part.clust <- subset_taxa(ps.genus, !(colnames(otu_table(ps.genus)) %in% colnames(otu_table(ps.genus.part.clust))))

#Add cluster abundance as one group
ps.genus.part.clust <- tax_glom(ps.genus.part.clust, taxrank = "Kingdom", NArm = F)
taxa_names(ps.genus.part.clust) <- "Partial_Cluster1"
ps.genus.sub.part.clust <- phyloseq(otu_table(ps.genus.sub.part.clust), tax_table(ps.genus.sub.part.clust), sample_data(ps.genus.sub.part.clust))
ps.genus.sub.part.clust <- merge_phyloseq(ps.genus.sub.part.clust, ps.genus.part.clust)
cat("\n","Phyloseq object with partial cluster 1 added:","\n")
ps.genus.sub.part.clust

#Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps.genus.sub.part.clust))
taxa.df <- rownames_to_column(taxa.df, "Sample.ID")

#Extract sample metadata from phyloseq object
meta.df <- data.frame(sample_data(ps.genus))
meta.df <- meta.df[,-1]
colnames(meta.df)[1] <- "Sample.ID"
cat("\n","Case/control included in analysis:","\n")
table(meta.df$case_control)

#Code fruits or vegetables daily variable
meta.df$fruits_or_vegetables_daily <- NA
meta.df$fruits_or_vegetables_daily[meta.df$fruits_or_vegetables == "At_least_once_a_day"] <- "Yes"
meta.df$fruits_or_vegetables_daily[meta.df$fruits_or_vegetables == "Few_times_a_month" | meta.df$fruits_or_vegetables == "Few_times_a_week"] <- "No"

#Grab some numbers to double-check Ns
cat("\n","sex:","\n")
table(meta.df$case_control, meta.df$sex)
cat("\n","Case age:","\n")
summary(subset(meta.df, case_control == "Case")$age)
cat("\n","Control age:","\n")
summary(subset(meta.df, case_control == "Control")$age)
cat("\n","location:","\n")
table(meta.df$case_control, meta.df$location)
cat("\n","Case stool travel time:","\n")
summary(subset(meta.df, case_control == "Case")$stool_travel_time)
cat("\n","Control stool travel time:","\n")
summary(subset(meta.df, case_control == "Control")$stool_travel_time)
cat("\n","digest problems day of stool collection:","\n")
table(meta.df$case_control, meta.df$ssdigest_prob)
cat("\n","fruits or vegetables daily:","\n")
table(meta.df$case_control, meta.df$fruits_or_vegetables_daily)

na.id <- meta.df[is.na(meta.df$sex) | is.na(meta.df$age) | is.na(meta.df$location) | is.na(meta.df$stool_travel_time) | is.na(meta.df$ssdigest_prob) | is.na(meta.df$fruits_or_vegetables_daily),]$Sample.ID
cat("\n","Total case/control in analysis after NA filtering:","\n")
table(meta.df[!(meta.df$Sample.ID %in% na.id),]$case_control)

#Set adjustment formula
adj.formula <- "sex+age+location+stool_travel_time+ssdigest_prob+fruits_or_vegetables_daily"

cat("\n","Model being tested:","\n")
paste("Genus ~ ","case_control+",adj.formula, sep="")

#Run ANCOM
ancom.results <- ANCOM.main(OTUdat=taxa.df,
                           Vardat=meta.df,
                           adjusted=TRUE,
                           repeated=FALSE,
                           main.var="case_control",
                           adj.formula=adj.formula,
                           repeat.var=NULL,
                           longitudinal=FALSE,
                           random.formula=NULL,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=1)

#add taxa designations to table
ancom.Wtaxa <- ancom.results$W.taxa
colnames(ancom.Wtaxa)[1] <- "Representative_ASV"
results <- dplyr::inner_join(ancom.Wtaxa, rownames_to_column(data.frame(tax_table(ps.genus.sub.part.clust)), "Representative_ASV"), by = "Representative_ASV")

write.table(results, "../Script_Output/Dataset1_Output/ANCOMv2_MWAS_Cluster1.txt", quote = F, sep="\t", row.names = F)

