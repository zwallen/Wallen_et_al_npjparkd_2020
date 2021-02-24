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
ps <- readRDS("../PhyloseqObjects/Dataset2/phyloseq.rds")

cat("\n","Phyloseq object summary:","\n")
ps

#Extract phylogenetic tree using the most abundant ASV as a root
phytree <- root(phy_tree(ps), names(sort(colSums(otu_table(ps)), decreasing=T))[1], resolve.root=T)
cat("\n", "Phylogenetic tree now rooted. Most abundant ASV used as root!","\n")

#Grab some numbers to double-check Ns
cat("\n", "Case/control included in analysis:","\n")
table(sample_data(ps)$case_control)
cat("\n","sex:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$sex, exclude=F)
cat("\n","Case age:","\n")
summary(subset(sample_data(ps), case_control == "Case")$age_GMQ)
cat("\n","Control age:","\n")
summary(subset(sample_data(ps), case_control == "Control")$age_GMQ)
cat("\n","Case bmi:","\n")
summary(subset(sample_data(ps), case_control == "Case")$bmi)
cat("\n","Control bmi:","\n")
summary(subset(sample_data(ps), case_control == "Control")$bmi)
cat("\n","loss 10lbs:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$loss_10lbs, exclude=F)
cat("\n","digest problems day of stool collection:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$ssdigest_prob, exclude=F)
cat("\n","p3m constipation:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$p3m_constipation, exclude=F)
cat("\n","alcohol currently:","\n")
table(sample_data(ps)$case_control, sample_data(ps)$do_you_drink_alcohol, exclude=F)
cat("\n","Case stool travel time:","\n")
summary(subset(sample_data(ps), case_control == "Case")$stool_travel_time)
cat("\n","Control stool travel time:","\n")
summary(subset(sample_data(ps), case_control == "Control")$stool_travel_time)

cat("\n","Model being tested:","\n")
paste("Distance ~ ","case_control+sex+age_GMQ+bmi+loss_10lbs+ssdigest_prob+p3m_constipation+do_you_drink_alcohol+stool_travel_time", sep="")

#filter out samples who are NAs for covariates
na.id <- sample_data(subset_samples(ps, is.na(sex) | is.na(age_GMQ) | is.na(bmi) | is.na(loss_10lbs) | is.na(ssdigest_prob) | is.na(p3m_constipation) | is.na(do_you_drink_alcohol) | is.na(stool_travel_time)))$HA_ID
ps <- subset_samples(ps, !(HA_ID %in% na.id))
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

#Run PERMANOVAs with Canberra, generalized unifrac
ps.canb <- adonis2(vegdist(otu_table(ps.ra), "canberra") ~ case_control+sex+age_GMQ+bmi+loss_10lbs+ssdigest_prob+p3m_constipation+do_you_drink_alcohol+stool_travel_time, data = as(sample_data(ps.ra), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(unifracs[, , "d_0.5"]) ~ case_control+sex+age_GMQ+bmi+loss_10lbs+ssdigest_prob+p3m_constipation+do_you_drink_alcohol+stool_travel_time, data = as(sample_data(ps.ra), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(vegdist(otu_table(ps.clr), "euclidean") ~ case_control+sex+age_GMQ+bmi+loss_10lbs+ssdigest_prob+p3m_constipation+do_you_drink_alcohol+stool_travel_time, data = as(sample_data(ps.clr), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

