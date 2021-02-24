#Compute SparCC correlations and emprical P-values for correlations
#Results used to generate Supplementary Figure 1, Dataset 1 cases
date()
library(phyloseq); packageVersion("phyloseq")
library(tibble); packageVersion("tibble")
library(dplyr); packageVersion("dplyr")

phyloseq.object <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(phyloseq.object, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps.genus))
taxa.df <- rownames_to_column(taxa.df, "SampleID")

#Extract sample metadata from phyloseq object
meta.df <- data.frame(sample_data(ps.genus))

#Merge taxa.df and meta.df
taxa.meta.df <- inner_join(taxa.df, meta.df, by = "SampleID")

#Format tables for SparCC and print so SparCC program can use them
case.table <- data.frame(t(subset(taxa.meta.df, case_control == "Case")[,2:dim(taxa.df)[2]]))
colnames(case.table) <- subset(taxa.meta.df, case_control == "Case")$SampleID
case.table <- inner_join(rownames_to_column(data.frame(tax_table(ps.genus)[,6]), "Taxa"), 
			 rownames_to_column(case.table, "Taxa"), by = "Taxa")[,-1]
colnames(case.table)[1] <- "OTU_id"

#Replace current labels with shorter taxonomy names
for (i in 1:length(taxa_names(ps.genus))){
        taxonomy <- as.vector(tax_table(ps.genus)[i, 1:6])
        taxonomy.sub <- taxonomy[!is.na(taxonomy)]
        if (length(taxonomy.sub) == 6){
                case.table[i,1] <- paste("g__",taxonomy.sub[6], sep="")
        }else if (length(taxonomy.sub) == 5){
                case.table[i,1] <- paste("f__",taxonomy.sub[5],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 4){
                case.table[i,1] <- paste("o__",taxonomy.sub[4],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 3){
                case.table[i,1] <- paste("c__",taxonomy.sub[3],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 2){
                case.table[i,1] <- paste("p__",taxonomy.sub[2],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 1){
                case.table[i,1] <- paste("k__",taxonomy.sub[1],"_unclass",sep="")
        }else{  
                case.table[i,1] <- "unclassified"
        }
}

write.table(case.table, "../Script_Output/Dataset1_Output/SparCC_case_table.txt", row.names=F, quote=F, sep="\t")

cat("\n","SparCC input data created","\n")

#Calculate SparCC correlation matrices
cat("\n","Beginning SparCC correlation calculations","\n")
system("SparCC.py ../Script_Output/Dataset1_Output/SparCC_case_table.txt --cor_file=../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Cases.txt --cov_file=../Script_Output/Dataset1_Output/SparCC_Cov_Matrix_Cases.txt")
cat("\n","Finished","\n")

#Create shuffled (with replacement) datasets for pseudo-pvalue calculation
cat("\n","Beginning construction of shuffled datasets","\n")
system("mkdir ../Script_Output/Dataset1_Output/SparCC_permutations_case")
system("MakeBootstraps.py ../Script_Output/Dataset1_Output/SparCC_case_table.txt -n 3000 --template=SparCC_case_table_perm#.txt --path=../Script_Output/Dataset1_Output/SparCC_permutations_case/")
cat("\n","Finished","\n")

#Calculate SparCC correlation matrices for each shuffled dataset
cat("\n","Beginning SparCC correlation calculations for shuffled datasets","\n")
system("cd ../Script_Output/Dataset1_Output/SparCC_permutations_case/; for i in SparCC_case_table_perm*.txt; do SparCC.py $i --cor_file=$i.cor.txt --cov_file=$i.cov.txt; done")
cat("\n","Finished","\n")

#Calculate pseudo-pvalues for SparCC correlations
cat("\n","Beginning pseudo-pvalue calculation for SparCC correlations","\n")
system("PseudoPvals.py ../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Cases.txt ../Script_Output/Dataset1_Output/SparCC_permutations_case/SparCC_case_table_perm#.txt.cor.txt 3000 -o ../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Cases_pvals.txt -t two_sided")
cat("\n","Finished","\n")

system("rm -r ../Script_Output/Dataset1_Output/SparCC_permutations_case")
system("rm ../Script_Output/Dataset1_Output/SparCC_case_table.txt")
