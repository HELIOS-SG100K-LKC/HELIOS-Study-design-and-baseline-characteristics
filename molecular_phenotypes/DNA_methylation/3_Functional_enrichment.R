##### Enrichment in regulatory regions

setwd("path/to/Epigenomics_Data/")

library(dplyr)
library(ggplot2)
library(cowplot)
library(ggnewscale)

##### Load result for each list and combine all to get results based on all 16,444 CpGs

####DHS
DHS <- read.delim("Results/CPG_list_00_DHS/CPG_list_00.850k.erc2-DHS.chart.tsv.gz")
DHS_Filenames <- unique(DHS$File)
DHS_random <- read.delim("Results/CPG_list_00_DHS/background.tsv.gz",header = F, col.names = c("File",paste0("Random",1:1000)))
DHS_random$File <- DHS_Filenames


List <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21")

for (i in List){
  Part_file <- read.delim(paste("Results/CPG_list_",i,"_DHS/CPG_list_",i,".850k.erc2-DHS.chart.tsv.gz",
                                sep = ""))
  DHS <- rbind(DHS,Part_file)
  Part_file_random <- read.delim(paste("Results/CPG_list_",i,"_DHS/background.tsv.gz",sep = ""),
                                 header = F, col.names = c("File",paste0("Random",1:1000)))
  Part_file_random$File <- DHS_Filenames
  DHS_random <- rbind(DHS_random,Part_file_random)
}

DHS_final <- as.data.frame(DHS %>% group_by(File,Accession,Cell,Tissue,Datatype) %>% summarise(
  Sentinel_Count = sum(TestStat),
  Random_Avg = sum(Random.set.Overlap)/1000))

DHS_random_final <- as.data.frame(DHS_random %>% group_by(File) %>%
                                    summarize(across(starts_with("Random"), ~sum(.x, na.rm = TRUE))))

rownames(DHS_random_final) <- DHS_random_final$File
DHS_random_final$File <- NULL
DHS_random_final_tr <- as.data.frame(t(DHS_random_final))

DHS_test <- DHS_final[,c(1,6)]
rownames(DHS_test) <- DHS_test$File
DHS_test$File <- NULL
DHS_test <- as.data.frame(t(DHS_test))

DHS_analytics <- c("File","Fold_Enrichment","Log2FE","Pval_Enrichment","Pval_Depletion","SD")
for (i in DHS_Filenames){
  Count <- as.numeric(DHS_test[[paste(i)]])
  Random_set <- as.vector(DHS_random_final_tr[[paste(i)]])
  Random_set_Avg <- mean(Random_set)
  SD_random <- sd(Random_set)
  enrichment_p <- sum(Random_set>=Count)/length(Random_set)
  depletion_p <- sum(Random_set<=Count)/length(Random_set)
  FE <- Count/Random_set_Avg
  log2FE <- log2(FE)
  Res <- c(paste(i),FE,log2FE,enrichment_p,depletion_p,SD_random)
  DHS_analytics <- rbind(DHS_analytics,Res)
}

DHS_analytics <- as.data.frame(DHS_analytics)
DHS_analytics <- DHS_analytics[-1,]
names(DHS_analytics) <- c("File","Fold_Enrichment","Log2FE","Pval_Enrichment","Pval_Depletion","SD")
DHS_analytics$Fold_Enrichment <- as.numeric(as.vector(DHS_analytics$Fold_Enrichment))
DHS_analytics$Log2FE <- as.numeric(as.vector(DHS_analytics$Log2FE))
DHS_analytics$Pval_Enrichment <- as.numeric(as.vector(DHS_analytics$Pval_Enrichment))
DHS_analytics$Pval_Depletion <- as.numeric(as.vector(DHS_analytics$Pval_Depletion))
rownames(DHS_analytics) <- NULL


DHS_final <- merge.data.frame(DHS_final,DHS_analytics,by="File")

Accessions <- DHS_final$Accession ## Needed to filter from Chromatin states


###### H3 ############

H3K4me1 <- read.delim("Results/CPG_list_00_H3K4me1/CPG_list_00.850k.erc2-H3K4me1.chart.tsv.gz",header = T)
H3K4me3 <- read.delim("Results/CPG_list_00_H3K4me3/CPG_list_00.850k.erc2-H3K4me3.chart.tsv.gz",header = T)
H3K9me3 <- read.delim("Results/CPG_list_00_H3K9me3/CPG_list_00.850k.erc2-H3K9me3.chart.tsv.gz",header = T)
H3K27me3 <- read.delim("Results/CPG_list_00_H3K27me3/CPG_list_00.850k.erc2-H3K27me3.chart.tsv.gz",header = T)
H3K36me3 <- read.delim("Results/CPG_list_00_H3K36me3/CPG_list_00.850k.erc2-H3K36me3.chart.tsv.gz",header = T)

H3_all <- rbind(H3K4me1,H3K4me3,H3K9me3,H3K27me3,H3K36me3)
H3_Filenames <- H3_all$File

H3K4me1_random <- read.delim("Results/CPG_list_00_H3K4me1/background.tsv.gz",header = F, col.names = c("File",paste0("Random",1:1000)))
H3K4me3_random <- read.delim("Results/CPG_list_00_H3K4me3/background.tsv.gz",header = F, col.names = c("File",paste0("Random",1:1000)))
H3K9me3_random <- read.delim("Results/CPG_list_00_H3K9me3/background.tsv.gz",header = F, col.names = c("File",paste0("Random",1:1000)))
H3K27me3_random <- read.delim("Results/CPG_list_00_H3K27me3/background.tsv.gz",header = F, col.names = c("File",paste0("Random",1:1000)))
H3K36me3_random <- read.delim("Results/CPG_list_00_H3K36me3/background.tsv.gz",header = F, col.names = c("File",paste0("Random",1:1000)))
H3_all_random <- rbind(H3K4me1_random,H3K4me3_random,H3K9me3_random,H3K27me3_random,H3K36me3_random)

H3_all_random$File <- H3_Filenames

for (i in List){
  Part_file_1 <- read.delim(paste("Results/CPG_list_",i,"_H3K4me1/CPG_list_",i,".850k.erc2-H3K4me1.chart.tsv.gz",sep = ""))
  Part_file_2 <- read.delim(paste("Results/CPG_list_",i,"_H3K4me3/CPG_list_",i,".850k.erc2-H3K4me3.chart.tsv.gz",sep = ""))
  Part_file_3 <- read.delim(paste("Results/CPG_list_",i,"_H3K9me3/CPG_list_",i,".850k.erc2-H3K9me3.chart.tsv.gz",sep = ""))
  Part_file_4 <- read.delim(paste("Results/CPG_list_",i,"_H3K27me3/CPG_list_",i,".850k.erc2-H3K27me3.chart.tsv.gz",sep = ""))
  Part_file_5 <- read.delim(paste("Results/CPG_list_",i,"_H3K36me3/CPG_list_",i,".850k.erc2-H3K36me3.chart.tsv.gz",sep = ""))
  Part_file <- rbind(Part_file_1,Part_file_2,Part_file_3,Part_file_4,Part_file_5)
  H3_all <- rbind(H3_all,Part_file)
  
  Part_file_random_1 <- read.delim(paste("Results/CPG_list_",i,"_H3K4me1/background.tsv.gz",sep = ""),
                                   header = F, col.names = c("File",paste0("Random",1:1000)))
  Part_file_random_2 <- read.delim(paste("Results/CPG_list_",i,"_H3K4me3/background.tsv.gz",sep = ""),
                                   header = F, col.names = c("File",paste0("Random",1:1000)))
  Part_file_random_3 <- read.delim(paste("Results/CPG_list_",i,"_H3K9me3/background.tsv.gz",sep = ""),
                                   header = F, col.names = c("File",paste0("Random",1:1000)))
  Part_file_random_4 <- read.delim(paste("Results/CPG_list_",i,"_H3K27me3/background.tsv.gz",sep = ""),
                                   header = F, col.names = c("File",paste0("Random",1:1000)))
  Part_file_random_5 <- read.delim(paste("Results/CPG_list_",i,"_H3K36me3/background.tsv.gz",sep = ""),
                                   header = F, col.names = c("File",paste0("Random",1:1000)))
  Part_file_random <- rbind(Part_file_random_1,Part_file_random_2,Part_file_random_3,Part_file_random_4,Part_file_random_5)
  Part_file_random$File <- H3_Filenames
  H3_all_random <- rbind(H3_all_random,Part_file_random)
}

H3_all_final <- as.data.frame(H3_all %>% group_by(File,Accession,Cell,Tissue,Datatype) %>% summarise(
  Sentinel_Count = sum(TestStat),
  Random_Avg = sum(Random.set.Overlap)/1000))

H3_all_random_final <- as.data.frame(H3_all_random %>% group_by(File) %>%
                                       summarize(across(starts_with("Random"), ~sum(.x, na.rm = TRUE))))

rownames(H3_all_random_final) <- H3_all_random_final$File
H3_all_random_final$File <- NULL
H3_all_random_final_tr <- as.data.frame(t(H3_all_random_final))

H3_all_test <- H3_all_final[,c(1,6)]
rownames(H3_all_test) <- H3_all_test$File
H3_all_test$File <- NULL
H3_all_test <- as.data.frame(t(H3_all_test))

H3_all_analytics <- c("File","Fold_Enrichment","Log2FE","Pval_Enrichment","Pval_Depletion","SD")
for (i in H3_Filenames){
  Count <- as.numeric(H3_all_test[[paste(i)]])
  Random_set <- as.vector(H3_all_random_final_tr[[paste(i)]])
  Random_set_Avg <- mean(Random_set)
  SD <- sd(Random_set)
  enrichment_p <- sum(Random_set>=Count)/length(Random_set)
  depletion_p <- sum(Random_set<=Count)/length(Random_set)
  FE <- Count/Random_set_Avg
  log2FE <- log2(FE)
  Res <- c(paste(i),FE,log2FE,enrichment_p,depletion_p,SD)
  H3_all_analytics <- rbind(H3_all_analytics,Res)
}

H3_all_analytics <- as.data.frame(H3_all_analytics)
H3_all_analytics <- H3_all_analytics[-1,]
names(H3_all_analytics) <- c("File","Fold_Enrichment","Log2FE","Pval_Enrichment","Pval_Depletion","SD")
H3_all_analytics$Fold_Enrichment <- as.numeric(as.vector(H3_all_analytics$Fold_Enrichment))
H3_all_analytics$Log2FE <- as.numeric(as.vector(H3_all_analytics$Log2FE))
H3_all_analytics$Pval_Enrichment <- as.numeric(as.vector(H3_all_analytics$Pval_Enrichment))
H3_all_analytics$Pval_Depletion <- as.numeric(as.vector(H3_all_analytics$Pval_Depletion))
rownames(H3_all_analytics) <- NULL


H3_all_final <- merge.data.frame(H3_all_final,H3_all_analytics,by="File")



######## Chromatin State ##############

CHR <- read.delim("Results/CPG_list_00_chromatin15state-all/CPG_list_00.850k.erc2-chromatin15state-all.chart.tsv.gz")
CHR_Filenames <- CHR$File
CHR_random <- read.delim("Results/CPG_list_00_chromatin15state-all/background.tsv.gz",header = F, col.names = c("File",paste0("Random",1:1000)))
CHR_random$File <- CHR_Filenames

List <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21")

for (i in List){
  Part_file <- read.delim(paste("Results/CPG_list_",i,"_chromatin15state-all/CPG_list_",i,".850k.erc2-chromatin15state-all.chart.tsv.gz",
                                sep = ""))
  CHR <- rbind(CHR,Part_file)
  Part_file_random <- read.delim(paste("Results/CPG_list_",i,"_chromatin15state-all/background.tsv.gz",sep = ""),
                                 header = F, col.names = c("File",paste0("Random",1:1000)))
  Part_file_random$File <- CHR_Filenames
  CHR_random <- rbind(CHR_random,Part_file_random)
}

CHR_final <- as.data.frame(CHR %>% group_by(File,Accession,Cell,Tissue,Datatype) %>% summarise(
  Sentinel_Count = sum(TestStat),
  Random_Avg = sum(Random.set.Overlap)/1000))

CHR_random_final <- as.data.frame(CHR_random %>% group_by(File) %>%
                                    summarize(across(starts_with("Random"), ~sum(.x, na.rm = TRUE))))

rownames(CHR_random_final) <- CHR_random_final$File
CHR_random_final$File <- NULL
CHR_random_final_tr <- as.data.frame(t(CHR_random_final))

CHR_test <- CHR_final[,c(1,6)]
rownames(CHR_test) <- CHR_test$File
CHR_test$File <- NULL
CHR_test <- as.data.frame(t(CHR_test))

CHR_analytics <- c("File","Fold_Enrichment","Log2FE","Pval_Enrichment","Pval_Depletion","SD")
for (i in CHR_Filenames){
  Count <- as.numeric(CHR_test[[paste(i)]])
  Random_set <- as.vector(CHR_random_final_tr[[paste(i)]])
  Random_set_Avg <- mean(Random_set)
  SD <- sd(Random_set)
  enrichment_p <- sum(Random_set>=Count)/length(Random_set)
  depletion_p <- sum(Random_set<=Count)/length(Random_set)
  FE <- Count/Random_set_Avg
  log2FE <- log2(FE)
  Res <- c(paste(i),FE,log2FE,enrichment_p,depletion_p,SD)
  CHR_analytics <- rbind(CHR_analytics,Res)
}

CHR_analytics <- as.data.frame(CHR_analytics)
CHR_analytics <- CHR_analytics[-1,]
names(CHR_analytics) <- c("File","Fold_Enrichment","Log2FE","Pval_Enrichment","Pval_Depletion","SD")
CHR_analytics$Fold_Enrichment <- as.numeric(as.vector(CHR_analytics$Fold_Enrichment))
CHR_analytics$Log2FE <- as.numeric(as.vector(CHR_analytics$Log2FE))
CHR_analytics$Pval_Enrichment <- as.numeric(as.vector(CHR_analytics$Pval_Enrichment))
CHR_analytics$Pval_Depletion <- as.numeric(as.vector(CHR_analytics$Pval_Depletion))
rownames(CHR_analytics) <- NULL



CHR_final <- merge.data.frame(CHR_final,CHR_analytics,by="File")
CHR_final <- CHR_final %>% filter(Accession %in% Accessions)


### Standardize Cell type names
Acc_Cell <- CHR_final[,c(2,3)]
names(Acc_Cell)[2] <- "Cell_Type"


All_Results <- rbind(DHS_final,H3_all_final,CHR_final)
All_Results <- merge.data.frame(Acc_Cell,All_Results,by="Accession")
All_Results <- All_Results[,c(3,1,5,2,6,7,8,9,10,11,12,13)]

#####Extract Final Results
write.table(All_Results,file = "Ehtnic_CpGs_All_Enrichment_results.txt",row.names = F,quote = F,sep = "\t")








################################################################################################################

 

#### TFBS Enrichment Analysis

TFBS_Count_file <- read.table("TF_Count_file.txt",col.names = c("TF","Background_Count","Sentinel_Count"))
TFBS_Count_file <- TFBS_Count_file %>% filter(Sentinel_Count>0)


Background_total <- 837222 # Total number of non-missing CpGs included in the analysis
Sentinel_total <- 16,444 ## Total number of ethnically diverse CpGs


TFBS_Count_file$FE <- (TFBS_Count_file$Sentinel_Count/Sentinel_total)/(TFBS_Count_file$Background_Count/Background_total)

TFBS_Count_file$Phyper <- phyper(TFBS_Count_file$Sentinel_Count-1,
                                  TFBS_Count_file$Background_Count,
                                  Background_total - TFBS_Count_file$Background_Count,
                                  Sentinel_total,
                                  lower.tail = F)

TFBS_Count_file$Padj <- p.adjust(TFBS_Count_file$phyper,method = "fdr")


write.table(TFBS_Count_file,file = TF_Enrichment_ethnic_CpGs.txt,row.names = F,quote = F,sep = " ")













