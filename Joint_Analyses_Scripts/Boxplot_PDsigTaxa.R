#Create boxplot of PD associated genera relative abundances
#Generates Figure 2
date()
library(phyloseq); packageVersion("phyloseq")
library(reshape2); packageVersion("reshape2")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(scales); packageVersion("scales")

target.taxa <- c("Porphyromonas", "Prevotella", "Corynebacterium_1", "Faecalibacterium", "Agathobacter", "Blautia", "Roseburia", "Butyricicoccus", "Fusicatenibacter", "Lachnospira", "Lachnospiraceae_ND3007_group", "Lachnospiraceae_UCG-004", "Oscillospira", "Bifidobacterium", "Lactobacillus")

### Dataset 1 ###
cat("\n", "Getting plot data ready for dataset 1","\n")
ps <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(ps, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Add count of 1 to asv table for plotting purposes
otu_table(ps.genus) <- otu_table(ps.genus)+1

#Transform counts to relative abundance
ps.genus.ra <- transform_sample_counts(ps.genus, function(x){x/sum(x)})
cat("\n","Genera counts transformed to proportions!","\n")

#Extract target taxa
ps.genus.filtered <- subset_taxa(ps.genus.ra, Genus %in% target.taxa)
cat("\n","Extracting target genera:","\n")
as.vector(tax_table(ps.genus.filtered)[,6])

#Get data ready for plotting
g.data <- data.frame(otu_table(ps.genus.filtered))
colnames(g.data) <- tax_table(ps.genus.filtered)[,6]
g.data <- g.data[target.taxa]
g.data <- cbind(data.frame(case_control=sample_data(ps.genus.filtered)$case_control), g.data, datasets=rep("Dataset 1", dim(g.data)[1]))
g.data.melt.1 <- melt(g.data)

### Dataset 2 ###
cat("\n", "Getting plot data ready for dataset 2","\n")
ps <- readRDS("../PhyloseqObjects/Dataset2/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(ps, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Add count of 1 to asv table for plotting purposes
otu_table(ps.genus) <- otu_table(ps.genus)+1

#Transform counts to relative abundance
ps.genus.ra <- transform_sample_counts(ps.genus, function(x){x/sum(x)})
cat("\n","Genera counts transformed to proportions!","\n")

#Extract target taxa
ps.genus.filtered <- subset_taxa(ps.genus.ra, Genus %in% target.taxa)
cat("\n","Extracting target genera:","\n")
as.vector(tax_table(ps.genus.filtered)[,6])

#Get data ready for plotting
g.data <- data.frame(otu_table(ps.genus.filtered))
colnames(g.data) <- tax_table(ps.genus.filtered)[,6]
g.data <- g.data[target.taxa]
g.data <- cbind(data.frame(case_control=sample_data(ps.genus.filtered)$case_control), g.data, datasets=rep("Dataset 2", dim(g.data)[1]))
g.data.melt.2 <- melt(g.data)

### Plotting ###

#Concatenate plot data for dataset1 and 2
g.data.melt <- rbind(g.data.melt.1, g.data.melt.2)

#Add extra data needed for plotting
g.data.melt$case_control_b <- recode(g.data.melt$case_control, Case = "a", Control = "b")

g.data.melt$cluster[g.data.melt$variable %in% target.taxa[1:3]] <- "Cluster 1"
g.data.melt$cluster[g.data.melt$variable %in% target.taxa[4:13]] <- "Cluster 2"
g.data.melt$cluster[g.data.melt$variable %in% target.taxa[14:15]] <- "Cluster 3"

#Create plot for relative abudances of tested genera
pdf("Boxplot_PDsigTaxa.pdf", height=10, width=20)
ggplot(data = g.data.melt, aes(x = variable, y = value, fill = case_control)) + 
  stat_boxplot(geom="errorbar", width=0.25, position=position_dodge(0.75)) +
  geom_boxplot(notch=T, outlier.size=-1) + 
  geom_dotplot(inherit.aes=F, aes(x = variable, y = value, color = case_control_b), binaxis = "y", stackdir = "center", binwidth = 0.03, position = position_dodge(0.75), dotsize = 1) +
  facet_grid(datasets~cluster, scale="free", space = "free", drop=TRUE) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 13)) + 
  theme(axis.text.y = element_text(size = 13), axis.title = element_text(size=17), legend.text = element_text(size = 13)) + 
  theme(strip.text = element_text(size = 20)) +
  labs(x="", y="Relative abundance (log10 scale)") + 
  scale_fill_manual(values=c("#00BFC4", "#E69F00")) + 
  scale_color_manual(values=c("#000000","#000000")) + 
  guides(fill=guide_legend(title=NULL), color=F)
dev.off()
system("mv Boxplot_PDsigTaxa.pdf ../Script_Output/Joint_Analyses_Output/")
