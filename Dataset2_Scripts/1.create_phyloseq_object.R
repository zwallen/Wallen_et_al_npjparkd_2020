###Create phyloseq object for dataset2
date()
library(phyloseq); packageVersion("phyloseq")

#Read in all data components
asv.table <- readRDS("../Script_Output/Dataset2_Output/seqtab.rds")
taxa.table <- readRDS("../Script_Output/Dataset2_Output/ASVTaxaSILVA132minBoot80.rds")
phylo.tree <- readRDS("../Script_Output/Dataset2_Output/seqtab_phylo_tree.rds")

metadata <- as.data.frame(readxl::read_xlsx("../Metadata/Dataset2_Metadata_6-09-20.xlsx", na = c("Unknown", "Not Sure"), .name_repair="minimal"))
rownames(metadata) <- metadata$HA_ID

#Create phyloseq object
phyloseq.object <- phyloseq(tax_table(taxa.table), sample_data(metadata),phy_tree(phylo.tree$tree),
                            otu_table(asv.table, taxa_are_rows = FALSE))
phyloseq.object

#Write phyloseq object to disk
saveRDS(phyloseq.object, "../PhyloseqObjects/Dataset2/phyloseq.rds")
