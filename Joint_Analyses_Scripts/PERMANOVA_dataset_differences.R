#Perform PERMANOVA for differences between datasets with and without adjustment for total sequence depth
#Results reported in manuscript text
#Note: results for PERMANOVA are printed straight to the terminal, so if running as job on SLURM scheduler results will be in the ".out" file in the "Output" directory
date()
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(microbiome); packageVersion("microbiome")

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

#Create combined phyloseq object, make sure all are collapsed to genus level
ps <- merge_phyloseq(dataset1.ps, dataset2.ps)
cat("\n","Joint phyloseq object summary:","\n")
ps

#Calculate total sequence depth
sample_data(ps)$total_seq_count <- rowSums(otu_table(ps))

#Get case control count
cat("\n","PD patient and control count for analysis:","\n")
table(sample_data(ps)$case_control)

#Get location count
sample_data(ps)$location[is.na(sample_data(ps)$location) & !is.na(sample_data(ps)$HA_ID)] <- "Birmingham,_AL"
cat("\n","Geographic location count for analysis:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$location)

#Create dataset variable
sample_data(ps)$datasets <- NA
sample_data(ps)$datasets[sample_data(ps)$location == "Birmingham,_AL"] <- "dataset2"
sample_data(ps)$datasets[sample_data(ps)$location != "Birmingham,_AL"] <- "dataset1"
cat("\n", "Dataset count for analysis:", "\n")
table(sample_data(ps)$case_control, sample_data(ps)$datasets)

#Transform counts with clr transformation
ps.clr <- transform(ps, transform="clr")
cat("\n","ASV counts transformed using centered log-ratio!","\n")

###Perform PERMANOVA for differences in datasets using Aitchison distance###

cat("\n","Running PERMANOVA using Aitchison distance...", "\n")

#Run PERMANOVA for dataset variable only
filt.ps <- subset_samples(ps.clr, !is.na(datasets))

datasets <- adonis2(vegdist(otu_table(filt.ps), "euclidean") ~ datasets, data = as(sample_data(filt.ps),"data.frame"), permutations = 99999)
cat("\n","PERMANOVA results for Dataset1 vs Dataset2:","\n")
datasets

#Run PERMANOVA for total seq count variable only
filt.ps <- subset_samples(ps.clr, !is.na(total_seq_count))

total_seq_count <- adonis2(vegdist(otu_table(filt.ps), "euclidean") ~ total_seq_count, data = as(sample_data(filt.ps),"data.frame"), permutations = 99999)
cat("\n","PERMANOVA results for total sequence count:","\n")
total_seq_count

#Run PERMANOVA for with dataset and total sequence count variables in model
filt.ps <- subset_samples(ps.clr, !is.na(datasets) & !is.na(total_seq_count))

datasets.total_seq_count <- adonis2(vegdist(otu_table(filt.ps), "euclidean") ~ datasets + total_seq_count, data = as(sample_data(filt.ps), "data.frame"), permutations = 99999, by = "margin")
cat("\n","PERMANOVA for Dataset1 vs Dataset2 and total sequence count:","\n")
datasets.total_seq_count

