library(RNOmni)
library(tidyverse)
library(data.table)


atc_codes <- paste0("^",c("A10B", "C10", "C0[2|3|7|8|9]"))
path_out <- "Med_pheno.txt"
longitudinal_path <- "detailed_longitudinal_data.gz"
pheno_path <- "covariate_file.txt.gz"
phenonames <- c("A10B", "C10", "BP_MED") #phenonames
min_age <- 10 #minimum aage requirement at the start of the follow-up (1995)
alive_only = "yes" #include only those who did not die before the start of the follow-up (1995)



  
  #read in the longitudinal data and phenofile
  print("Start drug phenotyper!!!")
  print("Reading in data...")
  print("Reading longitudinal data")
  dat <- fread(cmd=paste0("zcat ",longitudinal_path), data.table = F)
  dat <- subset(dat, SOURCE == "PURCH")
  print("Reading covariate data")
  pheno <- fread(cmd=paste0("zcat ",pheno_path), data.table = F)
  print("Reading done!")
  #use this if you want to limit by age, here alive and >= years old at the start of the drug-registy follow up, later you can subset to pheno <- subset(pheno, raw = 1)
  
  pheno_raw <- pheno[, c("FINNGENID","BL_YEAR","BL_AGE", "DEATH_AGE")]
  pheno_raw$AGE_1995 <- 1995 - (pheno_raw$BL_YEAR - pheno_raw$BL_AGE)
  pheno_raw$FU_TIME <- pheno_raw$DEATH_AGE-pheno_raw$AGE_1995
  if(min_age>0){
    pheno_raw <- subset(pheno_raw, AGE_1995>=min_age)
  }
  if(alive_only == "yes"){
    pheno_raw <- subset(pheno_raw, FU_TIME>0)
  }
  pheno_raw$raw <- 1
  pheno_raw <- pheno_raw[,c("FINNGENID","FU_TIME")]
  
  #pheno_raw <- merge(pheno, pheno_raw, by = "FINNGENID", all.x = T)
  #pheno <- pheno[,c("FINNGENID")]
  pheno <- left_join(pheno, pheno_raw, by ="FINNGENID")
  #here the atc-codes, use any and as many you want, start with "^" before the code to make sure to make an exact match
  
  print("Start endpointter!")
  
  for(i in 1:length(atc_codes)){
    print(i)
    print(paste0("Mining: ", atc_codes[i]))
    drug <- subset(dat, grepl(atc_codes[i], CODE1))
    first_first <- drug %>% arrange(EVENT_AGE)
    first_first <- first_first[1,c("FINNGENIED", "EVENT_AGE")]
    colnames(first_first)<-c("FINNGENID", "FIRST_AGE")
    last_first <- drug %>% arrange(EVENT_AGE)
    last_first <- last_first[1,c("FINNGENIED", "EVENT_AGE")]
    colnames(last_first)<- c("FINNGENID", "LAST_AGE")
    counts <- as.data.frame(drug %>% count(FINNGENID))
    print(nrow(counts))
    colnames(counts)[2]<- "drug_n"
    
    counts$drug_nUO <- counts$drug_n
    counts$drug_nUO <- RankNorm(counts$drug_nUO)
    
    drug <- left_join(pheno_raw, counts, by = "FINNGENID")
    drug <- left_join(drug, first_first, by = "FINNGENID")
    drug <- left_join(drug, last_first, by = "FINNGENID")
    drug$drug_n[is.na(drug$drug_n)]<-0
    drug$drug_count <- drug$drug_n
    drug$drug_n <- RankNorm(drug$drug_n)
    
    
    drug <- drug[,c("FINNGENID", "drug_n", "drug_nUO", "FIRST_AGE", "LAST_AGE")]
    colnames(drug)[2]<-phenonames[i]
    colnames(drug)[3]<-paste0(phenonames[i], "_UO")
    colnames(drug)[4]<-paste0(phenonames[i], "_count")
    colnames(drug)[5]<-paste0(phenonames[i], "_FIRST_AGE")
    colnames(drug)[6]<-paste0(phenonames[i], "_LAST_AGE")
    pheno <- merge(pheno, drug, by = "FINNGENID", all.x = T)
  }
  print("Endpointter finished!")
 
  
  
  #write out
  print("WRITING OUT!")
  write.table(pheno, path_out, quote=F, row.names=F, col.names=T, sep ="\t")
  system_cmd <- paste0("gzip ", path_out)
  system(system_cmd)

