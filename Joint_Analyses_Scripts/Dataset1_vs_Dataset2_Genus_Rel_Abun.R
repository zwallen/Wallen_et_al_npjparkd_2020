#Script that calculates relative abundance summary statistics for each genus in each dataset and puts them together into
#one convenient table

date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")

#Read in dataset1 phyloseq object
dataset1.ps <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")
cat("\n","dataset1 phyloseq object summary:","\n")
dataset1.ps

#Read in dataset2 phyloseq object
dataset2.ps <- readRDS("../PhyloseqObjects/Dataset2/phyloseq.rds")
cat("\n","dataset2 phyloseq object summary:","\n")
dataset2.ps

#Collapse to genus level
dataset1.genus <- tax_glom(dataset1.ps, taxrank="Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus for dataset1:","\n")
dataset1.genus

dataset2.genus <- tax_glom(dataset2.ps, taxrank="Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus for dataset2:","\n")
dataset2.genus

#Transform counts to relative abundances
dataset1.genus.ra <- transform_sample_counts(dataset1.genus, function(x) {x/sum(x)})

dataset2.genus.ra <- transform_sample_counts(dataset2.genus, function(x) {x/sum(x)})

###Loop through each genus and calculate average relative abundance for cases, controls, and the ratio for dataset1###
#Prep data for loop
ps <- dataset1.genus.ra
taxa.df <- data.frame(otu_table(ps))
taxa.df <- rownames_to_column(taxa.df,"SampleID")
meta.df <- data.frame(sample_data(ps))
cat("\n","dataset1 Case/Control included for analysis:","\n")
table(meta.df$case_control)
taxa.meta.df <- dplyr::inner_join(taxa.df, meta.df, by = "SampleID")

#Loop through genera calculating average case and control relative abundance and case/control ratio
dataset1.results <- data.frame()
for (i in 2:length(taxa.df[1,])){
  avg.case.prop <- mean(subset(taxa.meta.df, case_control == "Case")[,i])
  avg.control.prop <- mean(subset(taxa.meta.df, case_control == "Control")[,i])
  abundance.ratio <- avg.case.prop/avg.control.prop
  total.n <- sum(taxa.meta.df[,i] != 0)
  OTU <- colnames(taxa.meta.df)[i]
  dataset1.results <- rbind(dataset1.results, data.frame(dataset1_N = total.n, dataset1_CaseAvgRA = avg.case.prop, 
                                       dataset1_ContAvgRA = avg.control.prop, 
                                       dataset1_Ratio = abundance.ratio))
}
#Add taxonomy information
dataset1.results <- cbind(dataset1.results, data.frame(Kingdom = tax_table(dataset1.genus.ra)[,"Kingdom"],
                                     Phylum = tax_table(dataset1.genus.ra)[,"Phylum"],
                                     Class = tax_table(dataset1.genus.ra)[,"Class"],
                                     Order = tax_table(dataset1.genus.ra)[,"Order"],
                                     Family = tax_table(dataset1.genus.ra)[,"Family"],
                                     Genus = tax_table(dataset1.genus.ra)[,"Genus"]))

###Loop through each genus and calculate average relative abundance for cases, controls, and the ratio for dataset2###
#Prep data for loop
ps <- dataset2.genus.ra
taxa.df <- data.frame(otu_table(ps))
taxa.df <- rownames_to_column(taxa.df,"SampleID")
meta.df <- data.frame(sample_data(ps))
meta.df <- meta.df[,-1]
colnames(meta.df)[1] <- "SampleID"
cat("\n","dataset2 Case/Control included for analysis:","\n")
table(meta.df$case_control)
taxa.meta.df <- dplyr::inner_join(taxa.df, meta.df, by = "SampleID")

#Loop through genera calculating average case and control relative abundance and case/control ratio
dataset2.results <- data.frame()
for (i in 2:length(taxa.df[1,])){
  avg.case.prop <- mean(subset(taxa.meta.df, case_control == "Case")[,i])
  avg.control.prop <- mean(subset(taxa.meta.df, case_control == "Control")[,i])
  abundance.ratio <- avg.case.prop/avg.control.prop
  total.n <- sum(taxa.meta.df[,i] != 0)
  OTU <- colnames(taxa.meta.df)[i]
  dataset2.results <- rbind(dataset2.results, data.frame(dataset2_N = total.n, dataset2_CaseAvgRA = avg.case.prop, 
                                                 dataset2_ContAvgRA = avg.control.prop, 
                                                 dataset2_Ratio = abundance.ratio))
}
#Add taxonomy information
dataset2.results <- cbind(dataset2.results, data.frame(Kingdom = tax_table(dataset2.genus.ra)[,"Kingdom"],
                                               Phylum = tax_table(dataset2.genus.ra)[,"Phylum"],
                                               Class = tax_table(dataset2.genus.ra)[,"Class"],
                                               Order = tax_table(dataset2.genus.ra)[,"Order"],
                                               Family = tax_table(dataset2.genus.ra)[,"Family"],
                                               Genus = tax_table(dataset2.genus.ra)[,"Genus"]))

###Merge results and write table###
dataset1.dataset2.results <- merge(dataset1.results, dataset2.results, all=T)

write.table(dataset1.dataset2.results, "../Script_Output/Joint_Analyses_Output/Dataset1_vs_Dataset2_Genus_Rel_Abun.txt",
            quote = F, row.names = F, sep = '\t')
