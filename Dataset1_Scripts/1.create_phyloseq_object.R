###Create phyloseq object for dataset1
date()
library(phyloseq); packageVersion("phyloseq")

#Read in all data components
asv.table <- readRDS("../Script_Output/Dataset1_Output/seqtab.rds")
taxa.table <- readRDS("../Script_Output/Dataset1_Output/ASVTaxaSILVA132minBoot80.rds")
phylo.tree <- readRDS("../Script_Output/Dataset1_Output/seqtab_phylo_tree.rds")

metadata <- as.data.frame(readxl::read_xlsx("../Metadata/Dataset1_Metadata_3-08-20.xlsx", na = c("NA","not_applicable","not_provided","Unknown","Don't_Know"), .name_repair="minimal"))
rownames(metadata) <- metadata$SampleID

#Create phyloseq object
phyloseq.object <- phyloseq(tax_table(taxa.table), sample_data(metadata),phy_tree(phylo.tree$tree),
                            otu_table(asv.table, taxa_are_rows = FALSE))

#Make all data for FP0016201 and GMWA-1090 NA
sample_data(phyloseq.object)["10122.FP0016201",3:dim(sample_data(phyloseq.object))[2]] <- NA
sample_data(phyloseq.object)["10122.GMWA.1090",3:dim(sample_data(phyloseq.object))[2]] <- NA

#Put back data for case_control variable
sample_data(phyloseq.object)["10122.FP0016201","case_control"] <- "Case"
sample_data(phyloseq.object)["10122.GMWA.1090","case_control"] <- "Case"

cat("\n", "Phyloseq object summary:","\n")
phyloseq.object

#Write phyloseq object to disk
saveRDS(phyloseq.object, "../PhyloseqObjects/Dataset1/phyloseq.rds")
