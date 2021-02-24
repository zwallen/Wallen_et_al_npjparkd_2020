#Detect Genera significantly different between PD patients and controls using ANCOM adjusting for covariates
#Results reported in Table 2 and Supplementary Table 3
date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")
source("../Support_Files/ANCOM_updated_code.R")

ps <- readRDS("../PhyloseqObjects/Dataset2/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(ps, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps.genus))
taxa.df <- rownames_to_column(taxa.df, "Sample.ID")

#Extract sample metadata from phyloseq object
meta.df <- data.frame(sample_data(ps.genus))
meta.df <- meta.df[,-1]
colnames(meta.df)[1] <- "Sample.ID"
cat("\n","Case/control included in analysis:","\n")
table(meta.df$case_control)

#Grab some numbers to double-check Ns
cat("\n","sex:","\n")
table(meta.df$case_control, meta.df$sex)
cat("\n","Case age:","\n")
summary(subset(meta.df, case_control == "Case")$age_GMQ)
cat("\n","Control age:","\n")
summary(subset(meta.df, case_control == "Control")$age_GMQ)
cat("\n","Case bmi:","\n")
summary(subset(meta.df, case_control == "Case")$bmi)
cat("\n","Control bmi:","\n")
summary(subset(meta.df, case_control == "Control")$bmi)
cat("\n","p3m constipation:","\n")
table(meta.df$case_control, meta.df$p3m_constipation)

na.id <- meta.df[is.na(meta.df$sex) | is.na(meta.df$age_GMQ) | is.na(meta.df$bmi) | is.na(meta.df$p3m_constipation),]$Sample.ID
cat("\n","Total case/control in analysis after NA filtering:","\n")
table(meta.df[!(meta.df$Sample.ID %in% na.id),]$case_control)

#Set adjustment formula
adj.formula <- "sex+age_GMQ+bmi+p3m_constipation"

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
results <- dplyr::inner_join(ancom.Wtaxa, rownames_to_column(data.frame(tax_table(ps.genus)), "Representative_ASV"), by = "Representative_ASV")

write.table(results, "../Script_Output/Dataset2_Output/ANCOMv2_MWAS.txt", quote = F, sep="\t", row.names = F)

