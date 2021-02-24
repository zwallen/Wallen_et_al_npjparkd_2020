##Construct phylogenetic tree
library(DECIPHER)
library(phangorn)

asv.table <- readRDS("../../Script_Output/Dataset2_Output/seqtab.rds")
seqs <- colnames(asv.table)

names(seqs) <- seqs # This propagates to the tip labels of the tree

#Perform multiple sequence alignment
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
cat('\n', 'Multiple Sequence Alignment completed', '\n')

#Convert to object recognized by phangorn package
phang.align <- as.phyDat(as(alignment, "matrix"), type="DNA")
cat('\n','MSA converted to phyDat object','\n')

#Compute distance matrix
dm <- dist.ml(phang.align)
cat('\n','Distance matrix computed','\n')

#Perform neighbor joining tree construction
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
cat('\n', 'Phylogenetic tree contstructed','\n')

#Output tree to disk and save in Newick format
saveRDS(fitGTR, "../../Script_Output/Dataset2_Output/seqtab_phylo_tree.rds")
write.tree(fitGTR$tree,"../../Script_Output/Dataset2_Output/seqtab_phylo_tree.tre")

cat('\n','Done','\n')
