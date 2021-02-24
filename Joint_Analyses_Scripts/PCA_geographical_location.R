#Peform PCA with centered-log ratio transformation and color by geographical location
#Generates Figure 1
date()
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion("microbiome")
library(ggplot2); packageVersion("ggplot2")

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

#Create combined phyloseq object
ps <- merge_phyloseq(dataset1.ps, dataset2.ps)
cat("\n","Joint phyloseq object summary:","\n")
ps

#Get case control count
cat("\n","PD patient and control count for analysis:","\n")
table(sample_data(ps)$case_control)

#Get location count
sample_data(ps)$location[!is.na(sample_data(ps)$HA_ID)] <- "Birmingham,_AL"
cat("\n","Geographic location count for analysis:","\n")
table(sample_data(ps)$location)

#Exclude NA samples FP0016201 and GMWA-1090 for plotting purposes
ps <- subset_samples(ps, (SampleID != "10122.FP0016201" & SampleID != "10122.GMWA.1090") | is.na(SampleID))

### Cases and controls ###

#Transform counts to centered log-ratios
ps.clr <- transform(ps, transform="clr")
cat("\n", "ASV counts transformed to centered log-ratio","\n")

#Perform PCA
ps.ord <- ordinate(ps.clr, method="RDA", formula= ~1)

#Plot PCA
pdf("PCA_geographical_location.pdf", height=8, width=8)
plot_ordination(ps.clr, ps.ord, type="samples", axes=c(1,2), color="location", shape=NULL, label=NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
	axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.title = element_text(size=0), legend.text = element_text(size=20)) +
  geom_point(size=3) +
  scale_color_manual(breaks=c("Albany,_NY","Atlanta,_GA","Seattle,_WA","Birmingham,_AL"),
		      		   labels=c("Albany, NY","Atlanta, GA","Seattle, WA","Birmingham, AL"),
					   values=c("#F8766D","#7CAE00","#00CCFF","#C77CFF"),
					   na.translate = F)
dev.off()
system("mv PCA_geographical_location.pdf ../Script_Output/Joint_Analyses_Output/")

### Cases only ###

#Subset cases only
ps.sub <- subset_samples(ps, case_control == "Case")

#Transform counts to centered log-ratios
ps.clr <- transform(ps.sub, transform="clr")
cat("\n", "ASV counts transformed to centered log-ratio","\n")

#Perform PCA
ps.ord <- ordinate(ps.clr, method="RDA", formula= ~1)

#Plot PCA
pdf("PCA_geographical_location_cases.pdf", height=8, width=8)
plot_ordination(ps.clr, ps.ord, type="samples", axes=c(1,2), color="location", shape=NULL, label=NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
	axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.title = element_text(size=0), legend.text = element_text(size=20)) +
  geom_point(size=3) +
  scale_color_manual(breaks=c("Albany,_NY","Atlanta,_GA","Seattle,_WA","Birmingham,_AL"),
		      		   labels=c("Albany, NY","Atlanta, GA","Seattle, WA","Birmingham, AL"),
					   values=c("#F8766D","#7CAE00","#00CCFF","#C77CFF"),
					   na.translate = F)
dev.off()
system("mv PCA_geographical_location_cases.pdf ../Script_Output/Joint_Analyses_Output/")

### Controls only ###

#Subset controls only
ps.sub <- subset_samples(ps, case_control == "Control")

#Transform counts to centered log-ratios
ps.clr <- transform(ps.sub, transform="clr")
cat("\n", "ASV counts transformed to centered log-ratio","\n")

#Perform PCA
ps.ord <- ordinate(ps.clr, method="RDA", formula= ~1)

#Plot PCA
pdf("PCA_geographical_location_controls.pdf", height=8, width=8)
plot_ordination(ps.clr, ps.ord, type="samples", axes=c(1,2), color="location", shape=NULL, label=NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
	axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.title = element_text(size=0), legend.text = element_text(size=20)) +
  geom_point(size=3) +
  scale_color_manual(breaks=c("Albany,_NY","Atlanta,_GA","Seattle,_WA","Birmingham,_AL"),
		      		   labels=c("Albany, NY","Atlanta, GA","Seattle, WA","Birmingham, AL"),
					   values=c("#F8766D","#7CAE00","#00CCFF","#C77CFF"),
					   na.translate = F)
dev.off()
system("mv PCA_geographical_location_controls.pdf ../Script_Output/Joint_Analyses_Output/")

