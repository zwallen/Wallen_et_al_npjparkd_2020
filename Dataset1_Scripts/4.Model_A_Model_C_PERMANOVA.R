#Run PERMANOVA for case vs control not on certain PD medications
#Includes results reported for Model A and Model C, Table 1
#Note: results for PERMANOVA are printed straight to the terminal, so if running as job on SLURM scheduler results will be in the ".out" file in the "Output" directory
date()
library(phyloseq); packageVersion("phyloseq")
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

#Transform counts to proportions to normalize by library size. To be used in Canberra and Gen UniFrac calculation.
ps.ra <- transform_sample_counts(ps, function(x) {x/sum(x)})
cat("\n","ASV counts transformed to proportions for Canberra and Generalized UniFrac distance calculation.","\n")

#Transform counts to centered log-ratios. To be used in Aitchison distance calculations.
ps.clr <- transform(ps, transform="clr")
cat("\n", "ASV counts transformed to centered log-ratio for Aitchison distance calculation.","\n")

#Calculate UniFrac distances
gen.unifrac <- GUniFrac(otu_table(ps.ra), phytree, alpha=0.5)$unifracs[, , "d_0.5"]

#Calculate canberra distances
canberra <- as.matrix(vegdist(otu_table(ps.ra), "canberra"))

#Calculate Aitchison distances
aitchison <- as.matrix(vegdist(otu_table(ps.clr), "euclidean"))

### PD vs control ###

cat("\n","PERMANOVA for PD vs control","\n")
table(sample_data(ps)$case_control)

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(as.dist(canberra) ~ case_control, data = as(sample_data(ps), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(gen.unifrac) ~ case_control, data = as(sample_data(ps), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(as.dist(aitchison) ~ case_control, data = as(sample_data(ps), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

### PD not on drugs vs control ###

cat("\n","PERMANOVA for PD not on drugs vs control","\n")

#Subset distance matrices
ps.sub <- subset_samples(ps, (carbidopa_levodopa == "N" &
                   dopamine_agonist == "N" &
                   comt_inhibitor == "N" &
                   mao_b_inhibitor == "N" &
                   anticholinergic == "N" &
                   amantadine == "N") |
                   case_control == "Control")
gen.unifrac.sub <- gen.unifrac[rownames(gen.unifrac) %in% sample_names(ps.sub), colnames(gen.unifrac) %in% sample_names(ps.sub)]
canberra.sub <- canberra[rownames(canberra) %in% sample_names(ps.sub), colnames(canberra) %in% sample_names(ps.sub)]
aitchison.sub <- aitchison[rownames(aitchison) %in% sample_names(ps.sub), colnames(aitchison) %in% sample_names(ps.sub)]

table(sample_data(ps.sub)$case_control)

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(as.dist(canberra.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(gen.unifrac.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(as.dist(aitchison.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

### PD not on Ldopa vs control ###

cat("\n","PERMANOVA for PD not on Ldopa vs control","\n")

#Subset distance matrices
ps.sub <- subset_samples(ps, carbidopa_levodopa == "N" | case_control == "Control")
gen.unifrac.sub <- gen.unifrac[rownames(gen.unifrac) %in% sample_names(ps.sub), colnames(gen.unifrac) %in% sample_names(ps.sub)]
canberra.sub <- canberra[rownames(canberra) %in% sample_names(ps.sub), colnames(canberra) %in% sample_names(ps.sub)]
aitchison.sub <- aitchison[rownames(aitchison) %in% sample_names(ps.sub), colnames(aitchison) %in% sample_names(ps.sub)]

table(sample_data(ps.sub)$case_control)

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(as.dist(canberra.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(gen.unifrac.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(as.dist(aitchison.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

### PD not on COMT vs control ###

cat("\n","PERMANOVA for PD not on COMT vs control","\n")

#Subset distance matrices
ps.sub <- subset_samples(ps, comt_inhibitor == "N" | case_control == "Control")
gen.unifrac.sub <- gen.unifrac[rownames(gen.unifrac) %in% sample_names(ps.sub), colnames(gen.unifrac) %in% sample_names(ps.sub)]
canberra.sub <- canberra[rownames(canberra) %in% sample_names(ps.sub), colnames(canberra) %in% sample_names(ps.sub)]
aitchison.sub <- aitchison[rownames(aitchison) %in% sample_names(ps.sub), colnames(aitchison) %in% sample_names(ps.sub)]

table(sample_data(ps.sub)$case_control)

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(as.dist(canberra.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(gen.unifrac.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(as.dist(aitchison.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

### PD not on Anticholinergics vs control ###

cat("\n","PERMANOVA for PD not on Anticholinergics vs control","\n")

#Subset distance matrices
ps.sub <- subset_samples(ps, anticholinergic == "N" | case_control == "Control")
gen.unifrac.sub <- gen.unifrac[rownames(gen.unifrac) %in% sample_names(ps.sub), colnames(gen.unifrac) %in% sample_names(ps.sub)]
canberra.sub <- canberra[rownames(canberra) %in% sample_names(ps.sub), colnames(canberra) %in% sample_names(ps.sub)]
aitchison.sub <- aitchison[rownames(aitchison) %in% sample_names(ps.sub), colnames(aitchison) %in% sample_names(ps.sub)]

table(sample_data(ps.sub)$case_control)

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(as.dist(canberra.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(gen.unifrac.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(as.dist(aitchison.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

### PD not on MAO-B inhibitors vs control ###

cat("\n","PERMANOVA for PD not on MAO-B inhibitors vs control","\n")

#Subset distance matrices
ps.sub <- subset_samples(ps, mao_b_inhibitor == "N" | case_control == "Control")
gen.unifrac.sub <- gen.unifrac[rownames(gen.unifrac) %in% sample_names(ps.sub), colnames(gen.unifrac) %in% sample_names(ps.sub)]
canberra.sub <- canberra[rownames(canberra) %in% sample_names(ps.sub), colnames(canberra) %in% sample_names(ps.sub)]
aitchison.sub <- aitchison[rownames(aitchison) %in% sample_names(ps.sub), colnames(aitchison) %in% sample_names(ps.sub)]

table(sample_data(ps.sub)$case_control)

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(as.dist(canberra.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(gen.unifrac.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(as.dist(aitchison.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

### PD not on dopamine agonists vs control ###

cat("\n","PERMANOVA for PD not on dopamine agonists vs control","\n")

#Subset distance matrices
ps.sub <- subset_samples(ps, dopamine_agonist == "N" | case_control == "Control")
gen.unifrac.sub <- gen.unifrac[rownames(gen.unifrac) %in% sample_names(ps.sub), colnames(gen.unifrac) %in% sample_names(ps.sub)]
canberra.sub <- canberra[rownames(canberra) %in% sample_names(ps.sub), colnames(canberra) %in% sample_names(ps.sub)]
aitchison.sub <- aitchison[rownames(aitchison) %in% sample_names(ps.sub), colnames(aitchison) %in% sample_names(ps.sub)]

table(sample_data(ps.sub)$case_control)

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(as.dist(canberra.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(gen.unifrac.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(as.dist(aitchison.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

### PD not on amantadine vs control ###

cat("\n","PERMANOVA for PD not on amantadine vs control","\n")

#Subset distance matrices
ps.sub <- subset_samples(ps, amantadine == "N" | case_control == "Control")
gen.unifrac.sub <- gen.unifrac[rownames(gen.unifrac) %in% sample_names(ps.sub), colnames(gen.unifrac) %in% sample_names(ps.sub)]
canberra.sub <- canberra[rownames(canberra) %in% sample_names(ps.sub), colnames(canberra) %in% sample_names(ps.sub)]
aitchison.sub <- aitchison[rownames(aitchison) %in% sample_names(ps.sub), colnames(aitchison) %in% sample_names(ps.sub)]

table(sample_data(ps.sub)$case_control)

#Run PERMANOVAs with Canberra, generalized unifrac, and Aitchison distances
ps.canb <- adonis2(as.dist(canberra.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Canberra distance:","\n")
ps.canb

ps.gu <- adonis2(as.dist(gen.unifrac.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Generalized UniFrac distance:","\n")
ps.gu

ps.aitch <- adonis2(as.dist(aitchison.sub) ~ case_control, data = as(sample_data(ps.sub), "data.frame"), permutations = 99999, by = "margin")
cat("\n","Results for PERMANOVA with Aitchison distance:","\n")
ps.aitch

