#Run PERMANOVA for case vs control adjusting for covariates significantly different between cases and controls
#Includes results reported for Model B, Table 1
#Note: results for PERMANOVA are printed straight to the terminal, so if running as job on SLURM scheduler results will be in the ".out" file in the "Output" directory
date()
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
library(microbiome); packageVersion("microbiome")
library(vegan); packageVersion("vegan")
library(GUniFrac); packageVersion("GUniFrac")

#Read in phyloseq object
ps <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")

cat("\n","Phyloseq object summary:","\n")
ps

#Extract phylogenetic tree using the most abundant ASV as a root
phytree <- root(phy_tree(ps), names(sort(colSums(otu_table(ps)), decreasing=T))[1], resolve.root=T)
cat("\n", "Phylogenetic tree now rooted. Most abundant ASV used as root!","\n")

#Code fruits or vegetables daily variable
sample_data(ps)$fruits_or_vegetables_daily <- NA
sample_data(ps)$fruits_or_vegetables_daily[sample_data(ps)$fruits_or_vegetables == "At_least_once_a_day"] <- "Yes"
sample_data(ps)$fruits_or_vegetables_daily[sample_data(ps)$fruits_or_vegetables == "Few_times_a_month" | sample_data(ps)$fruits_or_vegetables == "Few_times_a_week"] <- "No"

#Code alcohol currently variable
sample_data(ps)$alcohol_currently <- NA
sample_data(ps)$alcohol_currently[sample_data(ps)$alcohol_amount == 0] <- "No"
sample_data(ps)$alcohol_currently[sample_data(ps)$alcohol_amount > 0] <- "Yes"

#Grab some numbers to double-check Ns
cat("\n", "Case/control included in analysis:","\n")
table(sample_data(ps)$case_control, exclude=F)
cat("\n","sex:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$sex, exclude=F)
cat("\n","Case age:","\n")
summary(subset(sample_data(ps), case_control == "Case")$age)
cat("\n","Control age:","\n")
summary(subset(sample_data(ps), case_control == "Control")$age)
cat("\n","Case bmi:","\n")
summary(subset(sample_data(ps), case_control == "Case")$bmi)
cat("\n","Control bmi:","\n")
summary(subset(sample_data(ps), case_control == "Control")$bmi)
cat("\n","location:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$location, exclude=F)
cat("\n","Case stool travel time:","\n")
summary(subset(sample_data(ps), case_control == "Case")$stool_travel_time)
cat("\n","Control stool travel time:","\n")
summary(subset(sample_data(ps), case_control == "Control")$stool_travel_time)
cat("\n","loss 10lbs:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$loss_10lbs, exclude=F)
cat("\n","digest problems day of stool collection:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$ssdigest_prob, exclude=F)
cat("\n","p3m constipation:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$p3m_constipation, exclude=F)
cat("\n","fruits or vegetables daily:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$fruits_or_vegetables_daily, exclude=F)
cat("\n","alcohol currently:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$alcohol_currently, exclude=F)

cat("\n","Model being tested:","\n")
paste("Distance ~ case_control+sex+age+bmi+location+stool_travel_time+loss_10lbs+ssdigest_prob+p3m_constipation+fruits_or_vegetables_daily+alcohol_currently")

#filter out samples who are NAs for covariates
na.id <- sample_data(subset_samples(ps, is.na(sex) | is.na(age) | is.na(bmi) | is.na(location) | is.na(stool_travel_time) | is.na(loss_10lbs) | is.na(ssdigest_prob) | is.na(p3m_constipation) | is.na(fruits_or_vegetables_daily) | is.na(alcohol_currently)))$SampleID
ps <- subset_samples(ps, !(SampleID %in% na.id))
cat("\n","Total case/control in analysis after NA filtering:","\n")
table(sample_data(ps)$case_control)

#Transform counts to proportions to normalize by library size. To be used in Canberra and Gen UniFrac calculation.
ps.ra <- transform_sample_counts(ps, function(x) {x/sum(x)})
cat("\n","ASV counts transformed to proportions for Canberra and Generalized UniFrac distance calculation.","\n")

#Transform counts to centered log-ratios. To be used in Aitchison distance calculations.
ps.clr <- transform(ps, transform="clr")
cat("\n", "ASV counts transformed to centered log-ratio for Aitchison distance calculation.","\n")

#Calculate UniFrac distances
unifracs <- GUniFrac(otu_table(ps.ra), phytree, alpha=0.5)$unifracs

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(vegdist(otu_table(ps.ra), "canberra") ~ case_control+sex+age+bmi+location+stool_travel_time+loss_10lbs+ssdigest_prob+p3m_constipation+fruits_or_vegetables_daily+alcohol_currently, data = as(sample_data(ps.ra), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(unifracs[, , "d_0.5"]) ~ case_control+sex+age+bmi+location+stool_travel_time+loss_10lbs+ssdigest_prob+p3m_constipation+fruits_or_vegetables_daily+alcohol_currently, data = as(sample_data(ps.ra), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(vegdist(otu_table(ps.clr), "euclidean") ~ case_control+sex+age+bmi+location+stool_travel_time+loss_10lbs+ssdigest_prob+p3m_constipation+fruits_or_vegetables_daily+alcohol_currently, data = as(sample_data(ps.clr), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

