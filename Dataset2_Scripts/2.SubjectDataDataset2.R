#Create summary statistics for metadata variables in Dataset 2 (Supplementary Table 1)
date()
library(phyloseq); packageVersion("phyloseq")

ps <- readRDS("../PhyloseqObjects/Dataset2/phyloseq.rds")

data <- data.frame(sample_data(ps))

#Create data frame for results
results <- data.frame()

#samples passing 16S QC
results <- rbind(results, data.frame(Variable = "Samples passing 16S QC", 
                                     PD_Case = table(data$case_control)[1],
                                     PD_Result = "-",
                                     Control = table(data$case_control)[2],
                                     Control_Result = "-",
                                     Total = length(na.omit(data$case_control)),
                                     Pvalue = "-"))

#Total unique ASVs detected
P.t <- length(sample_names(filter_taxa(subset_samples(ps, case_control == "Case"), function(x) sum(x) > 0, TRUE)))
P.tot <- length(taxa_names(filter_taxa(subset_samples(ps, case_control == "Case"), function(x) sum(x) > 0, TRUE)))
C.t <- length(sample_names(filter_taxa(subset_samples(ps, case_control == "Control"), function(x) sum(x) > 0, TRUE)))
C.tot <- length(taxa_names(filter_taxa(subset_samples(ps, case_control == "Control"), function(x) sum(x) > 0, TRUE)))

results <- rbind(results, data.frame(Variable = "Total Unique ASVs:", 
                                     PD_Case = P.t,
                                     PD_Result = P.tot,
                                     Control = C.t,
                                     Control_Result = C.tot,
                                     Total = P.t+C.t,
                                     Pvalue = "-"))

#Total genera detected
ps.genus <- tax_glom(ps, taxrank = "Genus", NArm = F)
P.t <- length(sample_names(filter_taxa(subset_samples(ps.genus, case_control == "Case"), function(x) sum(x) > 0, TRUE)))
P.tot <- length(taxa_names(filter_taxa(subset_samples(ps.genus, case_control == "Case"), function(x) sum(x) > 0, TRUE)))
C.t <- length(sample_names(filter_taxa(subset_samples(ps.genus, case_control == "Control"), function(x) sum(x) > 0, TRUE)))
C.tot <- length(taxa_names(filter_taxa(subset_samples(ps.genus, case_control == "Control"), function(x) sum(x) > 0, TRUE)))

results <- rbind(results, data.frame(Variable = "Total Genera:", 
                                     PD_Case = P.t,
                                     PD_Result = P.tot,
                                     Control = C.t,
                                     Control_Result = C.tot,
                                     Total = P.t+C.t,
                                     Pvalue = "-"))

#samples passing 16S and metadata QC
results <- rbind(results, data.frame(Variable = "Samples passing 16S and metadata QC", 
                                     PD_Case = table(data[!is.na(data$sex),]$case_control)[1],
                                     PD_Result = "-",
                                     Control = table(data[!is.na(data$sex),]$case_control)[2],
                                     Control_Result = "-",
                                     Total = length(na.omit(data$sex)),
                                     Pvalue = "-"))

#age
P.t <- length(na.omit(subset(data, case_control == "Case")$age_GMQ))
P.avg <- mean(na.omit(subset(data, case_control == "Case")$age_GMQ))
P.sd <- sd(na.omit(subset(data, case_control == "Case")$age_GMQ))
C.t <- length(na.omit(subset(data, case_control == "Control")$age_GMQ))
C.avg <- mean(na.omit(subset(data, case_control == "Control")$age_GMQ))
C.sd <- sd(na.omit(subset(data, case_control == "Control")$age_GMQ))
age.p <- wilcox.test(subset(data, case_control == "Case")$age_GMQ,
                     subset(data, case_control == "Control")$age_GMQ)$p.value
results <- rbind(results, data.frame(Variable = "Age: Mean;SD", 
                                     PD_Case = P.t,
                                     PD_Result = paste(round(P.avg, 1),round(P.sd, 1), sep = ";"),
                                     Control = C.t,
                                     Control_Result = paste(round(C.avg, 1),round(C.sd, 1), sep = ";"),
                                     Total = P.t+C.t,
                                     Pvalue = formatC(age.p, format = "e", digits = 1)))
#sex
P.f <- table(subset(data, case_control == "Case")$sex)[1]
P.m <- table(subset(data, case_control == "Case")$sex)[2]
C.f <- table(subset(data, case_control == "Control")$sex)[1]
C.m <- table(subset(data, case_control == "Control")$sex)[2]
sex.p <- fisher.test(matrix(c(P.f,P.m,C.f,C.m), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Sex: Male (% Male)", 
                                     PD_Case = P.f+P.m,
                                     PD_Result = paste(P.m," ","(",round(P.m/(P.f+P.m)*100, 0),"%",")", sep = ""),
                                     Control = C.f+C.m,
                                     Control_Result = paste(C.m," ","(",round(C.m/(C.f+C.m)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$sex)),
                                     Pvalue = formatC(sex.p, format = "e", digits = 1)))
#stool travel time
P.t <- length(na.omit(subset(data, case_control == "Case")$stool_travel_time))
P.avg <- mean(na.omit(subset(data, case_control == "Case")$stool_travel_time))
P.sd <- sd(na.omit(subset(data, case_control == "Case")$stool_travel_time))
C.t <- length(na.omit(subset(data, case_control == "Control")$stool_travel_time))
C.avg <- mean(na.omit(subset(data, case_control == "Control")$stool_travel_time))
C.sd <- sd(na.omit(subset(data, case_control == "Control")$stool_travel_time))
stool_travel_time.p <- wilcox.test(subset(data, case_control == "Case")$stool_travel_time,
                                   subset(data, case_control == "Control")$stool_travel_time)$p.value
results <- rbind(results, data.frame(Variable = "Stool sample travel time: Mean;SD", 
                                     PD_Case = P.t,
                                     PD_Result = paste(round(P.avg, 1),round(P.sd, 1), sep = ";"),
                                     Control = C.t,
                                     Control_Result = paste(round(C.avg, 1),round(C.sd, 1), sep = ";"),
                                     Total = P.t+C.t,
                                     Pvalue = formatC(stool_travel_time.p, format = "e", digits = 1)))
#race
P.n <- table(subset(data, case_control == "Case")$race)[1]
P.y <- table(subset(data, case_control == "Case")$race)[2]
C.n <- table(subset(data, case_control == "Control")$race)[1]
C.y <- table(subset(data, case_control == "Control")$race)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
race.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Race: White (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$race)),
                                     Pvalue = formatC(race.p, format = "e", digits = 1)))

#BMI
P.t <- length(na.omit(subset(data, case_control == "Case")$bmi))
P.avg <- mean(na.omit(subset(data, case_control == "Case")$bmi))
P.sd <- sd(na.omit(subset(data, case_control == "Case")$bmi))
C.t <- length(na.omit(subset(data, case_control == "Control")$bmi))
C.avg <- mean(na.omit(subset(data, case_control == "Control")$bmi))
C.sd <- sd(na.omit(subset(data, case_control == "Control")$bmi))
bmi.p <- wilcox.test(subset(data, case_control == "Case")$bmi,
                     subset(data, case_control == "Control")$bmi)$p.value
results <- rbind(results, data.frame(Variable = "BMI: Mean;SD", 
                                     PD_Case = P.t,
                                     PD_Result = paste(round(P.avg, 1),round(P.sd, 1), sep = ";"),
                                     Control = C.t,
                                     Control_Result = paste(round(C.avg, 1),round(C.sd, 1), sep = ";"),
                                     Total = P.t+C.t,
                                     Pvalue = formatC(bmi.p, format = "e", digits = 1)))
#lost > 10lbs last year
P.n <- table(subset(data, case_control == "Case")$loss_10lbs)[1]
P.y <- table(subset(data, case_control == "Case")$loss_10lbs)[2]
C.n <- table(subset(data, case_control == "Control")$loss_10lbs)[1]
C.y <- table(subset(data, case_control == "Control")$loss_10lbs)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
loss_10lbs.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Lost 10 pounds in Past Year: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$loss_10lbs)),
                                     Pvalue = formatC(loss_10lbs.p, format = "e", digits = 1)))
#gained > 10lbs last year
P.n <- table(subset(data, case_control == "Case")$gained_10lbs)[1]
P.y <- table(subset(data, case_control == "Case")$gained_10lbs)[2]
C.n <- table(subset(data, case_control == "Control")$gained_10lbs)[1]
C.y <- table(subset(data, case_control == "Control")$gained_10lbs)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
gained_10lbs.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Gained 10 pounds in Past Year: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$gained_10lbs)),
                                     Pvalue = formatC(gained_10lbs.p, format = "e", digits = 1)))

#do you drink alcohol
P.n <- table(subset(data, case_control == "Case")$do_you_drink_alcohol)[1]
P.y <- table(subset(data, case_control == "Case")$do_you_drink_alcohol)[2]
C.n <- table(subset(data, case_control == "Control")$do_you_drink_alcohol)[1]
C.y <- table(subset(data, case_control == "Control")$do_you_drink_alcohol)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
do_you_drink_alcohol.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Drinks alcohol: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$do_you_drink_alcohol)),
                                     Pvalue = formatC(do_you_drink_alcohol.p, format = "e", digits = 1)))

#do you smoke
P.n <- table(subset(data, case_control == "Case")$do_you_smoke)[1]
P.y <- table(subset(data, case_control == "Case")$do_you_smoke)[2]
C.n <- table(subset(data, case_control == "Control")$do_you_smoke)[1]
C.y <- table(subset(data, case_control == "Control")$do_you_smoke)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
do_you_smoke.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Tobacco: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$do_you_smoke)),
                                     Pvalue = formatC(do_you_smoke.p, format = "e", digits = 1)))
#caffeine
P.n <- table(subset(data, case_control == "Case")$do_you_drink_caffeinated_beverages)[1]
P.y <- table(subset(data, case_control == "Case")$do_you_drink_caffeinated_beverages)[2]
C.n <- table(subset(data, case_control == "Control")$do_you_drink_caffeinated_beverages)[1]
C.y <- table(subset(data, case_control == "Control")$do_you_drink_caffeinated_beverages)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
do_you_drink_caffeinated_beverages.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Caffeine: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$do_you_drink_caffeinated_beverages)),
                                     Pvalue = formatC(do_you_drink_caffeinated_beverages.p, format = "e", digits = 1)))

#No bowel movement for 3 days or more prior to stool collection
P.n <- table(subset(data, case_control == "Case")$ss_constipation)[1]
P.y <- table(subset(data, case_control == "Case")$ss_constipation)[2]
C.n <- table(subset(data, case_control == "Control")$ss_constipation)[1]
C.y <- table(subset(data, case_control == "Control")$ss_constipation)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
ss_constipation.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Constipation 3 Days or More Prior to Stool Collection: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$ss_constipation)),
                                     Pvalue = formatC(ss_constipation.p, format = "e", digits = 1)))

#diarrhea day of stool collection
P.n <- table(subset(data, case_control == "Case")$ss_diarrhea)[1]
P.y <- table(subset(data, case_control == "Case")$ss_diarrhea)[2]
C.n <- table(subset(data, case_control == "Control")$ss_diarrhea)[1]
C.y <- table(subset(data, case_control == "Control")$ss_diarrhea)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
ss_diarrhea.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Diarrhea Day of Stool Collection: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$ss_diarrhea)),
                                     Pvalue = formatC(ss_diarrhea.p, format = "e", digits = 1)))

#abdominal pain day of stool collection
P.n <- table(subset(data, case_control == "Case")$ssabdominal_pain)[1]
P.y <- table(subset(data, case_control == "Case")$ssabdominal_pain)[2]
C.n <- table(subset(data, case_control == "Control")$ssabdominal_pain)[1]
C.y <- table(subset(data, case_control == "Control")$ssabdominal_pain)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
ssabdominal_pain.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "GI Pain Day of Stool Collection: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$ssabdominal_pain)),
                                     Pvalue = formatC(ssabdominal_pain.p, format = "e", digits = 1)))

#excess gas day of stool collection
P.n <- table(subset(data, case_control == "Case")$ss_excess_gas)[1]
P.y <- table(subset(data, case_control == "Case")$ss_excess_gas)[2]
C.n <- table(subset(data, case_control == "Control")$ss_excess_gas)[1]
C.y <- table(subset(data, case_control == "Control")$ss_excess_gas)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
ss_excess_gas.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Excess Gas Day of Stool Collection: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$ss_excess_gas)),
                                     Pvalue = formatC(ss_excess_gas.p, format = "e", digits = 1)))
#bloating day of stool collection
P.n <- table(subset(data, case_control == "Case")$ss_bloating)[1]
P.y <- table(subset(data, case_control == "Case")$ss_bloating)[2]
C.n <- table(subset(data, case_control == "Control")$ss_bloating)[1]
C.y <- table(subset(data, case_control == "Control")$ss_bloating)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
ss_bloating.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Bloating Day of Stool Collection: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$ss_bloating)),
                                     Pvalue = formatC(ss_bloating.p, format = "e", digits = 1)))

#Gastrointestinal discomfort on day of collection
P.n <- table(subset(data, case_control == "Case")$ssdigest_prob)[1]
P.y <- table(subset(data, case_control == "Case")$ssdigest_prob)[2]
C.n <- table(subset(data, case_control == "Control")$ssdigest_prob)[1]
C.y <- table(subset(data, case_control == "Control")$ssdigest_prob)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
ssdigest_prob.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "GI Discomfort Day Of Stool Collection: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$ssdigest_prob)),
                                     Pvalue = formatC(ssdigest_prob.p, format = "e", digits = 1)))

#constipation past 3 months
P.n <- table(subset(data, case_control == "Case")$p3m_constipation)[1]
P.y <- table(subset(data, case_control == "Case")$p3m_constipation)[2]
C.n <- table(subset(data, case_control == "Control")$p3m_constipation)[1]
C.y <- table(subset(data, case_control == "Control")$p3m_constipation)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
p3m_constipation.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Constipation in past 3 months: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$p3m_constipation)),
                                     Pvalue = formatC(p3m_constipation.p, format = "e", digits = 1)))
#diarrhea past 3 months
P.n <- table(subset(data, case_control == "Case")$p3m_diarrhea)[1]
P.y <- table(subset(data, case_control == "Case")$p3m_diarrhea)[2]
C.n <- table(subset(data, case_control == "Control")$p3m_diarrhea)[1]
C.y <- table(subset(data, case_control == "Control")$p3m_diarrhea)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
p3m_diarrhea.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Diarrhea in past 3 months: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$p3m_diarrhea)),
                                     Pvalue = formatC(p3m_diarrhea.p, format = "e", digits = 1)))

#colitis
P.n <- table(subset(data, case_control == "Case")$colitis)[1]
P.y <- table(subset(data, case_control == "Case")$colitis)[2]
C.n <- table(subset(data, case_control == "Control")$colitis)[1]
C.y <- table(subset(data, case_control == "Control")$colitis)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
colitis.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Colitis: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$colitis)),
                                     Pvalue = formatC(colitis.p, format = "e", digits = 1)))

#IBS
P.n <- table(subset(data, case_control == "Case")$irrit_bowel_syndrome)[1]
P.y <- table(subset(data, case_control == "Case")$irrit_bowel_syndrome)[2]
C.n <- table(subset(data, case_control == "Control")$irrit_bowel_syndrome)[1]
C.y <- table(subset(data, case_control == "Control")$irrit_bowel_syndrome)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
irrit_bowel_syndrome.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "IBS: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$irrit_bowel_syndrome)),
                                     Pvalue = formatC(irrit_bowel_syndrome.p, format = "e", digits = 1)))

#Crohn's disease
P.n <- table(subset(data, case_control == "Case")$crohns_disease)[1]
P.y <- table(subset(data, case_control == "Case")$crohns_disease)[2]
C.n <- table(subset(data, case_control == "Control")$crohns_disease)[1]
C.y <- table(subset(data, case_control == "Control")$crohns_disease)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
crohns_disease.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Crohn's disease: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$crohns_disease)),
                                     Pvalue = formatC(crohns_disease.p, format = "e", digits = 1)))
#IBD
P.n <- table(subset(data, case_control == "Case")$inflam_bowel_disease)[1]
P.y <- table(subset(data, case_control == "Case")$inflam_bowel_disease)[2]
C.n <- table(subset(data, case_control == "Control")$inflam_bowel_disease)[1]
C.y <- table(subset(data, case_control == "Control")$inflam_bowel_disease)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
inflam_bowel_disease.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "IBD: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$inflam_bowel_disease)),
                                     Pvalue = formatC(inflam_bowel_disease.p, format = "e", digits = 1)))
#Ulcers
P.n <- table(subset(data, case_control == "Case")$ulcers)[1]
P.y <- table(subset(data, case_control == "Case")$ulcers)[2]
C.n <- table(subset(data, case_control == "Control")$ulcers)[1]
C.y <- table(subset(data, case_control == "Control")$ulcers)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
ulcers.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Ulcers: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$ulcers)),
                                     Pvalue = formatC(ulcers.p, format = "e", digits = 1)))
#SIBO
P.n <- table(subset(data, case_control == "Case")$SIBO)[1]
P.y <- table(subset(data, case_control == "Case")$SIBO)[2]
C.n <- table(subset(data, case_control == "Control")$SIBO)[1]
C.y <- table(subset(data, case_control == "Control")$SIBO)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
SIBO.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "SIBO: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$SIBO)),
                                     Pvalue = formatC(SIBO.p, format = "e", digits = 1)))
#celiac disease
P.n <- table(subset(data, case_control == "Case")$celiac_disease)[1]
P.y <- table(subset(data, case_control == "Case")$celiac_disease)[2]
C.n <- table(subset(data, case_control == "Control")$celiac_disease)[1]
C.y <- table(subset(data, case_control == "Control")$celiac_disease)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
celiac_disease.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Celiac disease: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$celiac_disease)),
                                     Pvalue = formatC(celiac_disease.p, format = "e", digits = 1)))
#cancer of digestive system
P.n <- table(subset(data, case_control == "Case")$GI_cancer)[1]
P.y <- table(subset(data, case_control == "Case")$GI_cancer)[2]
C.n <- table(subset(data, case_control == "Control")$GI_cancer)[1]
C.y <- table(subset(data, case_control == "Control")$GI_cancer)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
GI_cancer.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "GI Cancer: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$GI_cancer)),
                                     Pvalue = formatC(GI_cancer.p, format = "e", digits = 1)))
#Any of above GI problems (crohns, IBD, IBS, colitis, ulcers, SIBO, celiac disease, GI cancer)
P.n <- table(subset(data, case_control == "Case")$intest_disease)[1]
P.y <- table(subset(data, case_control == "Case")$intest_disease)[2]
C.n <- table(subset(data, case_control == "Control")$intest_disease)[1]
C.y <- table(subset(data, case_control == "Control")$intest_disease)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
intest_disease.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Intestinal Disease: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$intest_disease)),
                                     Pvalue = formatC(intest_disease.p, format = "e", digits = 1)))

#currently taking antibiotics
P.n <- table(subset(data, case_control == "Case")$antibiotics_bool)[1]
P.y <- table(subset(data, case_control == "Case")$antibiotics_bool)[2]
C.n <- table(subset(data, case_control == "Control")$antibiotics_bool)[1]
C.y <- table(subset(data, case_control == "Control")$antibiotics_bool)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
antibiotics_bool.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Currently taking antibiotics: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$antibiotics_bool)),
                                     Pvalue = formatC(antibiotics_bool.p, format = "e", digits = 1)))
#Taken antibiotics in past 3 months
P.n <- table(subset(data, case_control == "Case")$p3m_Antibiotics_bool)[1]
P.y <- table(subset(data, case_control == "Case")$p3m_Antibiotics_bool)[2]
C.n <- table(subset(data, case_control == "Control")$p3m_Antibiotics_bool)[1]
C.y <- table(subset(data, case_control == "Control")$p3m_Antibiotics_bool)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
if (is.na(C.n)){C.n <- 0}
if (is.na(C.y)){C.y <- 0}
p3m_Antibiotics_bool.p <- fisher.test(matrix(c(P.n,P.y,C.n,C.y), nrow = 2))$p.value
results <- rbind(results, data.frame(Variable = "Taken antibiotics in past 3 months: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = C.n+C.y,
                                     Control_Result = paste(C.y," ","(",round(C.y/(C.n+C.y)*100, 0),"%",")", sep = ""),
                                     Total = length(na.omit(data$p3m_Antibiotics_bool)),
                                     Pvalue = formatC(p3m_Antibiotics_bool.p, format = "e", digits = 1)))

#Disease duration
P.t <- length(na.omit(subset(data, case_control == "Case")$parkinson_disease_duration_GMQ))
P.avg <- mean(na.omit(subset(data, case_control == "Case")$parkinson_disease_duration_GMQ))
P.sd <- sd(na.omit(subset(data, case_control == "Case")$parkinson_disease_duration_GMQ))
results <- rbind(results, data.frame(Variable = "Disease duration: Mean;SD", 
                                     PD_Case = P.t,
                                     PD_Result = paste(round(P.avg, 1),round(P.sd, 1), sep = ";"),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = P.t,
                                     Pvalue = "-"))

#Carbidopa/levodopa
P.n <- table(subset(data, case_control == "Case")$levodopa)[1]
P.y <- table(subset(data, case_control == "Case")$levodopa)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
results <- rbind(results, data.frame(Variable = "Carbidopa/Levodopa: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = length(na.omit(data$levodopa)),
                                     Pvalue = "-"))
#Levodopa dose
P.t <- length(na.omit(subset(data, case_control == "Case")$total_levodopa_dose))
P.avg <- mean(na.omit(subset(data, case_control == "Case")$total_levodopa_dose))
P.sd <- sd(na.omit(subset(data, case_control == "Case")$total_levodopa_dose))
results <- rbind(results, data.frame(Variable = "Levodopa Dose: Mean;SD", 
                                     PD_Case = P.t,
                                     PD_Result = paste(round(P.avg, 1),round(P.sd, 1), sep = ";"),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = P.t,
                                     Pvalue = "-"))
#Dopamine agonist
P.n <- table(subset(data, case_control == "Case")$dopamine_agonist)[1]
P.y <- table(subset(data, case_control == "Case")$dopamine_agonist)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
results <- rbind(results, data.frame(Variable = "Dopamine agonist: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = length(na.omit(data$dopamine_agonist)),
                                     Pvalue = "-"))
#MAO-B inhibitor
P.n <- table(subset(data, case_control == "Case")$mao_b_inhibitor)[1]
P.y <- table(subset(data, case_control == "Case")$mao_b_inhibitor)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
results <- rbind(results, data.frame(Variable = "MAO-B inhibitor: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = length(na.omit(data$mao_b_inhibitor)),
                                     Pvalue = "-"))
#Amantadine
P.n <- table(subset(data, case_control == "Case")$amantadine)[1]
P.y <- table(subset(data, case_control == "Case")$amantadine)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
results <- rbind(results, data.frame(Variable = "Amantadine: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = length(na.omit(data$amantadine)),
                                     Pvalue = "-"))

#COMT inhibitor
P.n <- table(subset(data, case_control == "Case")$comt_inhibitor)[1]
P.y <- table(subset(data, case_control == "Case")$comt_inhibitor)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
results <- rbind(results, data.frame(Variable = "COMT inhibitor: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = length(na.omit(data$comt_inhibitor)),
                                     Pvalue = "-"))

#Anticholinergics
P.n <- table(subset(data, case_control == "Case")$anticholinergic)[1]
P.y <- table(subset(data, case_control == "Case")$anticholinergic)[2]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
results <- rbind(results, data.frame(Variable = "Anticholinergics: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = length(na.omit(data$anticholinergic)),
                                     Pvalue = "-"))
#Unmedicated
P.n <- table(subset(data, case_control == "Case")$pd_medications)[2]
P.y <- table(subset(data, case_control == "Case")$pd_medications)[1]
if (is.na(P.n)){P.n <- 0}
if (is.na(P.y)){P.y <- 0}
results <- rbind(results, data.frame(Variable = "Unmedicated: N (%)", 
                                     PD_Case = P.n+P.y,
                                     PD_Result = paste(P.y," ","(",round(P.y/(P.n+P.y)*100, 0),"%",")", sep = ""),
                                     Control = "-",
                                     Control_Result = "-",
                                     Total = length(na.omit(data$pd_medications)),
                                     Pvalue = "-"))

write.table(results, "../Script_Output/Dataset2_Output/SubjectDataDataset2.txt",
            row.names = F, quote = F, sep = '\t')


