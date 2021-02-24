#Create node and edge file from SparCC output to use in Gephi
#Used as input to Gephi to create Supplementary Figure 1, Dataset 1 cases and controls
date()
library(phyloseq); packageVersion("phyloseq")

phyloseq.object <- readRDS("../PhyloseqObjects/Dataset1/phyloseq.rds")

#Collapse taxa to Genus level
ps.genus <- tax_glom(phyloseq.object, taxrank = "Genus", NArm = F)
cat("\n","Collapsing ASVs by Genus:", "\n")
ps.genus

#Transform counts to relative abundance
ps.genus.ra <- transform_sample_counts(ps.genus, function(x){x/sum(x)})
cat("\n","Genera counts transformed to proportions!","\n")

ps.case <- subset_samples(ps.genus.ra, case_control == "Case")
ps.cont <- subset_samples(ps.genus.ra, case_control == "Control")

case.cor <- read.table("../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Cases.txt", header = T, row.names = "OTU_id", stringsAsFactors = F)
cont.cor <- read.table("../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Controls.txt", header = T, row.names = "OTU_id", stringsAsFactors = F)
case.p <- read.table("../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Cases_pvals.txt", header = T, row.names = "OTU_id", stringsAsFactors = F)
cont.p <- read.table("../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Controls_pvals.txt", header = T, row.names = "OTU_id", stringsAsFactors = F)

#Calculate mean relative abundance for calculating mean relative abundance ratios
avg.case.prop <- c()
avg.control.prop <- c()
for (i in 1:dim(otu_table(ps.genus.ra))[2]){
    
  avg.case.prop <- append(avg.case.prop, mean(otu_table(ps.case)[,i]))
  avg.control.prop <- append(avg.control.prop, mean(otu_table(ps.cont)[,i]))
 
}

#Create node file (only need one because case and control node files should be identical)
nodes.case <- data.frame(Id=rownames(case.cor), 
                    Label=rownames(case.cor),
                    Kingdom=as.vector(tax_table(ps.case)[,1]),
                    Phylum=as.vector(tax_table(ps.case)[,2]),
                    Class=as.vector(tax_table(ps.case)[,3]),
                    Order=as.vector(tax_table(ps.case)[,4]),
                    Family=as.vector(tax_table(ps.case)[,5]),
                    Genus=as.vector(tax_table(ps.case)[,6]),
		    FC=avg.case.prop/avg.control.prop, 
                    row.names = NULL)
nodes.cont <- data.frame(Id=rownames(cont.cor), 
                    Label=rownames(cont.cor),
                    Kingdom=as.vector(tax_table(ps.cont)[,1]),
                    Phylum=as.vector(tax_table(ps.cont)[,2]),
                    Class=as.vector(tax_table(ps.cont)[,3]),
                    Order=as.vector(tax_table(ps.cont)[,4]),
                    Family=as.vector(tax_table(ps.cont)[,5]),
                    Genus=as.vector(tax_table(ps.cont)[,6]),
		    FC=avg.case.prop/avg.control.prop, 
                    row.names = NULL)
if (identical(nodes.case, nodes.cont)){
	nodes <- nodes.case
}else{
	stop("Node files for cases and controls do not match. They should be identical!")
}

#Create edge files
corr.direction <- c()
corr.direction[case.cor[lower.tri(case.cor)] > 0] <- "+"
corr.direction[case.cor[lower.tri(case.cor)] < 0] <- "-"
edges.case <- data.frame(Source=t(combn(rownames(case.cor), 2))[,1],
                         Target=t(combn(rownames(case.cor), 2))[,2],
                         Weight=abs(case.cor[lower.tri(case.cor)]),
			 CorrDirection=corr.direction,
                         SparCCpval=case.p[lower.tri(case.p)]
			 )
corr.direction <- c()
corr.direction[cont.cor[lower.tri(cont.cor)] > 0] <- "+"
corr.direction[cont.cor[lower.tri(cont.cor)] < 0] <- "-"
edges.cont <- data.frame(Source=t(combn(rownames(cont.cor), 2))[,1],
                         Target=t(combn(rownames(cont.cor), 2))[,2],
                         Weight=abs(cont.cor[lower.tri(cont.cor)]),
			 CorrDirection=corr.direction,
                         SparCCpval=cont.p[lower.tri(cont.p)]
			 )

#Get list of PD associated genera
ancom.results.dataset1 <- read.table("../Script_Output/Dataset1_Output/ANCOMv2_MWAS.txt", header=T)
ancom.results.dataset2 <- read.table("../Script_Output/Dataset2_Output/ANCOMv2_MWAS.txt", header=T)
ancom.results <- dplyr::inner_join(ancom.results.dataset1, ancom.results.dataset2, by=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxa.ids <- c()
for (i in 1:dim(ancom.results)[1]){
        taxonomy <- as.vector(ancom.results[i, c("Kingdom","Phylum","Class","Order","Family","Genus")])
        taxonomy.sub <- taxonomy[!is.na(taxonomy)]
        if (length(taxonomy.sub) == 6){
                taxa.id <- paste("g__",taxonomy.sub[6], sep="")
        }else if (length(taxonomy.sub) == 5){
                taxa.id <- paste("f__",taxonomy.sub[5],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 4){
                taxa.id <- paste("o__",taxonomy.sub[4],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 3){
                taxa.id <- paste("c__",taxonomy.sub[3],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 2){
                taxa.id <- paste("p__",taxonomy.sub[2],"_unclass",sep="")
        }else if (length(taxonomy.sub) == 1){
                taxa.id <- paste("k__",taxonomy.sub[1],"_unclass",sep="")
        }else{
                taxa.id <- "unclassified"
        }
	taxa.ids <- append(taxa.ids, taxa.id)
}
ancom.results$taxa.ids <- taxa.ids
target.taxa <- ancom.results[ancom.results$detected_0.8.x == "TRUE" & ancom.results$detected_0.8.y == "TRUE",]$taxa.ids

#Tag nodes corresponding to PD associated genera denoting if they are increased or decreased in PD
nodes$PD_sig[nodes$Id %in% target.taxa & nodes$FC > 1] <- "Yes_Increased"
nodes$PD_sig[nodes$Id %in% target.taxa & nodes$FC < 1] <- "Yes_Decreased"
nodes$PD_sig[is.na(nodes$PD_sig)] <- "No"

#Tag edges containing PD associated genera
edges.case$PD_sig[edges.case$Source %in% target.taxa | edges.case$Target %in% target.taxa] <- "Yes"
edges.case$PD_sig[is.na(edges.case$PD_sig)] <- "No"
edges.cont$PD_sig[edges.cont$Source %in% target.taxa | edges.cont$Target %in% target.taxa] <- "Yes"
edges.cont$PD_sig[is.na(edges.cont$PD_sig)] <- "No"

#write out files
write.csv(nodes, "../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_nodeTable.csv",
          row.names = F, quote = F)
write.csv(edges.case, "../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Cases_edgeTable.csv",
          row.names = F, quote = F)
write.csv(edges.cont, "../Script_Output/Dataset1_Output/SparCC_Corr_Matrix_Controls_edgeTable.csv",
          row.names = F, quote = F)

