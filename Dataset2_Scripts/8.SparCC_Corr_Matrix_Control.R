#Compute SparCC correlations and emprical P-values for correlations
#Results used to generate Supplementary Figure 1, Dataset 2 controls
date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")
library(dplyr); packageVersion("dplyr")

phyloseq.object <- readRDS("../PhyloseqObjects/Dataset2/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(phyloseq.object, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps.genus))
taxa.df <- rownames_to_column(taxa.df, "SampleID")

#Extract sample metadata from phyloseq object
meta.df <- data.frame(sample_data(ps.genus))
meta.df <- meta.df[,-1]
colnames(meta.df)[1] <- "SampleID"

#Merge taxa.df and meta.df
taxa.meta.df <- inner_join(taxa.df, meta.df, by = "SampleID")

#Format tables for SparCC and print so SparCC program can use them
cont.table <- data.frame(t(subset(taxa.meta.df, case_control == "Control")[,2:dim(taxa.df)[2]]))
colnames(cont.table) <- subset(taxa.meta.df, case_control == "Control")$SampleID
cont.table <- inner_join(rownames_to_column(data.frame(tax_table(ps.genus)[,6]), "Taxa"),
			 rownames_to_column(cont.table, "Taxa"), by = "Taxa")[,-1]
colnames(cont.table)[1] <- "OTU_id"

#Replace current labels with shorter taxonomy names
for (i in 1:length(taxa_names(ps.genus))){
        taxonomy <- as.vector(tax_table(ps.genus)[i, 1:6])
        taxonomy.sub <- taxonomy[!is.na(taxonomy)]
        if (length(taxonomy.sub) == 6){
                cont.table[i,1] <- paste("g__",taxonomy.sub[6], sep="")
        }else if (length(taxonomy.sub) == 5){
                cont.table[i,1] <- paste("f__",taxonomy.sub[5],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 4){
                cont.table[i,1] <- paste("o__",taxonomy.sub[4],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 3){
                cont.table[i,1] <- paste("c__",taxonomy.sub[3],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 2){
                cont.table[i,1] <- paste("p__",taxonomy.sub[2],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 1){
                cont.table[i,1] <- paste("k__",taxonomy.sub[1],"_unclass",sep="")
        }else{
                cont.table[i,1] <- "unclassified"
        }
}

write.table(cont.table, "../Script_Output/Dataset2_Output/SparCC_cont_table.txt", row.names=F, quote=F, sep="\t")

cat("\n","SparCC input data created","\n")

#Calculate SparCC correlation matrices
cat("\n","Beginning SparCC correlation calculations","\n")
system("SparCC.py ../Script_Output/Dataset2_Output/SparCC_cont_table.txt --cor_file=../Script_Output/Dataset2_Output/SparCC_Corr_Matrix_Controls.txt --cov_file=../Script_Output/Dataset2_Output/SparCC_Cov_Matrix_Controls.txt")
cat("\n","Finished","\n")

#Create shuffled (with replacement) datasets for pseudo-pvalue calculation (WARNING: THIS REQUIRES A LOT OF SPACE!!!)
cat("\n","Beginning construction of shuffled datasets","\n")
system("mkdir ../Script_Output/Dataset2_Output/SparCC_permutations_cont")
system("MakeBootstraps.py ../Script_Output/Dataset2_Output/SparCC_cont_table.txt -n 3000 --template=SparCC_cont_table_perm#.txt --path=../Script_Output/Dataset2_Output/SparCC_permutations_cont/")
cat("\n","Finished","\n")

#Calculate SparCC correlation matrices for each shuffled dataset
cat("\n","Beginning SparCC correlation calculations for shuffled datasets","\n")
system("cd ../Script_Output/Dataset2_Output/SparCC_permutations_cont/; for i in SparCC_cont_table_perm*.txt; do SparCC.py $i --cor_file=$i.cor.txt --cov_file=$i.cov.txt; done")
cat("\n","Finished","\n")

#Calculate pseudo-pvalues for SparCC correlations
cat("\n","Beginning pseudo-pvalue calculation for SparCC correlations","\n")
system("PseudoPvals.py ../Script_Output/Dataset2_Output/SparCC_Corr_Matrix_Controls.txt ../Script_Output/Dataset2_Output/SparCC_permutations_cont/SparCC_cont_table_perm#.txt.cor.txt 3000 -o ../Script_Output/Dataset2_Output/SparCC_Corr_Matrix_Controls_pvals.txt -t two_sided")
cat("\n","Finished","\n")

system("rm -r ../Script_Output/Dataset2_Output/SparCC_permutations_cont")
system("rm ../Script_Output/Dataset2_Output/SparCC_cont_table.txt")

