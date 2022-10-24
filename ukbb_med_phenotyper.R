rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(RNOmni)




ukb <- fread("zcat  ukbb_covariate_file.tsv.gz") #load the UKBB covariatefile that includes the whole dataset
colnames(ukb)[1]<-"eid"

dat <- fread("zcat ukbb_gp_prescriptions_file.txt.gz") #load the gp scripts file
colnames(dat)[1]<-"eid"

dat$bnf_code <- gsub('\\.', '', dat$bnf_code) #take dots out of bnf codes

allcounts <- as.data.frame(dat%>% count(eid)) #count the number of individuals included in drug purchase data


bfs <- c('^0000', '^000[0|1|2]|^0001','^000002') #bnf codes for 3 different phenotypes 1), 2), 3), have to be specified
phenonames <- c("phenoname1", "phenoname2", "phenoname3") #corresponding atc codes for 1), 2), 3), have to be specified


pheno <- ukb

#loop through the bnf codes to get the counts for 1, 2, 3 on all individuals, include zeros or not =_UO extension and normalize
#and merge with phenofile

for(i in 1:length(bfs)){
  print(i)
  drug <- subset(dat, grepl(bfs[i], bnf_code))
  counts <- as.data.frame(drug %>% count(eid))
  print(nrow(counts))
  colnames(counts)[2]<- "drug_n"
  drug <- merge(allcounts, counts, by = "eid", all.x = T)
  drug$drug_nUO <- drug$drug_n
  drug$drug_n[is.na(drug$drug_n)]<-0
  drug <- drug[,c("eid", "drug_n", "drug_nUO")]
  colnames(drug)[2]<-phenonames[i]
  colnames(drug)[3]<-paste0(phenonames[i], "_UO")
  pheno <- merge(pheno, drug, by = "eid", all.x = T)
}

out_name <- "X" #specify X f.e. med_file.txt

write.table(pheno, "UKB_drug_counts_pheno.txt", quote=F, row.names=F, col.names=T, sep ="\t")