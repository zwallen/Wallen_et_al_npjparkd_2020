#Merge individual run seqence variant tables, and remove chimeras
library(dada2); packageVersion("dada2")
library(dplyr); packageVersion("dplyr")

#merge ASV tables from multiple runs
table1 <- readRDS("../../Script_Output/Dataset2_Output/seqtab01.rds")
table2 <- readRDS("../../Script_Output/Dataset2_Output/seqtab02.rds")
table3 <- readRDS("../../Script_Output/Dataset2_Output/seqtab03.rds")
table4 <- readRDS("../../Script_Output/Dataset2_Output/seqtab04.rds")
table5 <- readRDS("../../Script_Output/Dataset2_Output/seqtab05.rds")
table6 <- readRDS("../../Script_Output/Dataset2_Output/seqtab06.rds")
table.all <- mergeSequenceTables(table1, table2, table3, table4, table5, table6)
cat("\n","Sequence variant tables merged.","\n")

#remove chimeras
seqtab <- removeBimeraDenovo(table.all, method="consensus", multithread=TRUE)
cat("\n","Fraction of chimeras:", sum(seqtab)/sum(table.all), "\n")

#write merged, chimera free table to disk
saveRDS(seqtab, "../../Script_Output/Dataset2_Output/seqtab.rds")

#write ASV table
seqtab.df <- add_rownames(as.data.frame(seqtab), "SampleID")
write.table(seqtab.df, "../../Script_Output/Dataset2_Output/seqtab.txt", sep="\t", row.names=F, quote=F)

