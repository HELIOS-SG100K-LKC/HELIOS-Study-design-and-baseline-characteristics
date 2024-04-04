setwd("~/Path/to/Genomic_Data")

## packages

library(dplyr)
library(ggplot2)
library(cowplot)

##### Load all Complete Phenotype and covariate files

Phenotype_All <- read.csv("Batch-normImpDataAll_missing25_zerovar_1073metabolites_235repVis2Excl_9pc1OutliersExcl.csv",header = T)


Covar_file <- read.table("HELIOS_Covar_file.txt",header = T) ### Age and Gender Information
names(Covar_file) <- c("IID","Age","Gender")
SSID_batch <- read.table("HELIOS_ssid.txt",header = T)
Covar_file <- merge.data.frame(Covar_file,SSID_batch,by="IID")


######## Inverse Normal Transformation Function (Preserves missing values as NA)

inverse_normal_transform <- function(x) {
  non_na_indices <- which(!is.na(x))
  n <- length(non_na_indices)
  a <- 3/8
  ranks <- rank(x[non_na_indices])
  transformed_values <- x
  transformed_values[non_na_indices] <- qnorm((ranks - a) / (n + 1 - 2*a))
  return(transformed_values)
}

########################



## get list of phenotypes
Pheno_list <- names(Phenotype_All)
Pheno_list <- Pheno_list[-1]


##  Prepare Ethnicty Specific Phenotype File
#Chinese
Chinese_IDs <- read.table("Chinese.ids.txt",header = T)
Chinese_Phenotype_all <- Phenotype_All %>% filter(IID %in% Chinese_IDs$IID)
Chinese_Phenotype_all$FID <- Chinese_Phenotype_all$IID
Chinese_Phenotype_all <- Chinese_Phenotype_all[,c(1075,1:1074)]




##### Remove outliers (5 sd) ; Followed by Inverse Normal Trasformation

for (i in Pheno_list){
  mean_chi <- mean(Chinese_Phenotype_all[[paste0(i)]],na.rm = T)
  sd_chi  <- sd(Chinese_Phenotype_all[[paste0(i)]],na.rm = T)
  
  Chinese_Phenotype_all[[paste0(i)]] <- ifelse(Chinese_Phenotype_all[[paste0(i)]] > mean_chi + (5*sd_chi) | 
                                                 Chinese_Phenotype_all[[paste0(i)]] < mean_chi - (5*sd_chi),
                                               NA,Chinese_Phenotype_all[[paste0(i)]])
  
  Chinese_Phenotype_all[[paste0(i)]] <- inverse_normal_transform(Chinese_Phenotype_all[[paste0(i)]])
}


### Distribution Plots
pdf("Chinese_Pheno_distribution_raw_transformed.pdf", height = 10, width = 20)
for (i in Pheno_list) {
  p_raw <- ggplot(Chinese_Phenotype_all_raw, aes(x = .data[[i]])) +
    geom_density(alpha = 0.5, fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 8, face = "bold"))
  
  p_outrem_trans <- ggplot(Chinese_Phenotype_all, aes(x = .data[[i]])) +
    geom_density(alpha = 0.5, fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 12, face = "bold"))
  
  combined_plot <- plot_grid(p_raw, p_outrem_trans, labels = c("Raw Data", "Outlier Removed & INT"), ncol = 2)
  
  print(combined_plot)
}
dev.off()


## Split data into multiple files for regenie [each file with FID, IID and 36 Phenotypes]. 
#This will help run Regenie on each file for ease of computation.

IDs <- Chinese_Phenotype_all[,c(1,2)]
chunks <- ceiling(1073/36) ## 1073 phenotypes & 36 phenotypes in each file

for (i in 1:chunks) {
  # Calculate the columns for the current chunk
  start_col <- 3 + (i - 1) * 36
  end_col <- min(2 + i * 36, 1075) ## 1075 is the total number of columns
  
  # Select columns for the current chunk
  current_chunk <- as.data.frame(cbind(IDs, Chinese_Phenotype_all[, c(start_col:end_col)]))
  
  # Save the current chunk to a file (adjust the filename as needed)
  write.table(current_chunk, file = paste0("Phenotypes/Chinese_Phenotype_File_", i, ".txt"), sep = " ", row.names = FALSE,quote = F)
}


####Covariate File for Chinese

Chinese_Covar <- Covar_file %>% filter(IID %in% Chinese_Phenotype_all$IID)
Chinese_PCs <- read.table("HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06.eigenvec",header=T)
Chinese_Covar <- merge.data.frame(Chinese_Covar,Chinese_PCs,by="IID")
Chinese_Covar <- Chinese_Covar[,c(5,1:4,6:24)] ## FID, IID, Age, Sex, Batch, PC1-20
write.table(Chinese_Covar,file="Chinese_metab_covar.txt",row.names = F,quote = F,sep = " ")




#######
##  Prepare Ethnicty Specific Phenotype File
#Indian
Indian_IDs <- read.table("Indian.ids.txt",header = T)
Indian_Phenotype_all <- Phenotype_All %>% filter(IID %in% Indian_IDs$IID)
Indian_Phenotype_all$FID <- Indian_Phenotype_all$IID
Indian_Phenotype_all <- Indian_Phenotype_all[,c(1075,1:1074)]




##### Remove outliers (5 sd) ; Followed by Inverse Normal Trasformation

for (i in Pheno_list){
  mean_chi <- mean(Indian_Phenotype_all[[paste0(i)]],na.rm = T)
  sd_chi  <- sd(Indian_Phenotype_all[[paste0(i)]],na.rm = T)
  
  Indian_Phenotype_all[[paste0(i)]] <- ifelse(Indian_Phenotype_all[[paste0(i)]] > mean_chi + (5*sd_chi) | 
                                                Indian_Phenotype_all[[paste0(i)]] < mean_chi - (5*sd_chi),
                                              NA,Indian_Phenotype_all[[paste0(i)]])
  
  Indian_Phenotype_all[[paste0(i)]] <- inverse_normal_transform(Indian_Phenotype_all[[paste0(i)]])
}


### Distribution Plots
pdf("Indian_Pheno_distribution_raw_transformed.pdf", height = 10, width = 20)
for (i in Pheno_list) {
  p_raw <- ggplot(Indian_Phenotype_all_raw, aes(x = .data[[i]])) +
    geom_density(alpha = 0.5, fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 8, face = "bold"))
  
  p_outrem_trans <- ggplot(Indian_Phenotype_all, aes(x = .data[[i]])) +
    geom_density(alpha = 0.5, fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 12, face = "bold"))
  
  combined_plot <- plot_grid(p_raw, p_outrem_trans, labels = c("Raw Data", "Outlier Removed & INT"), ncol = 2)
  
  print(combined_plot)
}
dev.off()


## Split data into multiple files for regenie [each file with FID, IID and 36 Phenotypes]. 
#This will help run Regenie on each file for ease of computation.

IDs <- Indian_Phenotype_all[,c(1,2)]
chunks <- ceiling(1073/36) ## 1073 phenotypes & 36 phenotypes in each file

for (i in 1:chunks) {
  # Calculate the columns for the current chunk
  start_col <- 3 + (i - 1) * 36
  end_col <- min(2 + i * 36, 1075) ## 1075 is the total number of columns
  
  # Select columns for the current chunk
  current_chunk <- as.data.frame(cbind(IDs, Indian_Phenotype_all[, c(start_col:end_col)]))
  
  # Save the current chunk to a file (adjust the filename as needed)
  write.table(current_chunk, file = paste0("Phenotypes/Indian_Phenotype_File_", i, ".txt"), sep = " ", row.names = FALSE,quote = F)
}


####Covariate File for Indian

Indian_Covar <- Covar_file %>% filter(IID %in% Indian_Phenotype_all$IID)
Indian_PCs <- read.table("HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06.eigenvec",header=T)
Indian_Covar <- merge.data.frame(Indian_Covar,Indian_PCs,by="IID")
Indian_Covar <- Indian_Covar[,c(5,1:4,6:24)] ## FID, IID, Age, Sex, Batch, PC1-20
write.table(Indian_Covar,file="Indian_metab_covar.txt",row.names = F,quote = F,sep = " ")



#######
##  Prepare Ethnicty Specific Phenotype File
#Malay
Malay_IDs <- read.table("Malay.ids.txt",header = T)
Malay_Phenotype_all <- Phenotype_All %>% filter(IID %in% Malay_IDs$IID)
Malay_Phenotype_all$FID <- Malay_Phenotype_all$IID
Malay_Phenotype_all <- Malay_Phenotype_all[,c(1075,1:1074)]




##### Remove outliers (5 sd) ; Followed by Inverse Normal Trasformation

for (i in Pheno_list){
  mean_chi <- mean(Malay_Phenotype_all[[paste0(i)]],na.rm = T)
  sd_chi  <- sd(Malay_Phenotype_all[[paste0(i)]],na.rm = T)
  
  Malay_Phenotype_all[[paste0(i)]] <- ifelse(Malay_Phenotype_all[[paste0(i)]] > mean_chi + (5*sd_chi) | 
                                               Malay_Phenotype_all[[paste0(i)]] < mean_chi - (5*sd_chi),
                                             NA,Malay_Phenotype_all[[paste0(i)]])
  
  Malay_Phenotype_all[[paste0(i)]] <- inverse_normal_transform(Malay_Phenotype_all[[paste0(i)]])
}


### Distribution Plots
pdf("Malay_Pheno_distribution_raw_transformed.pdf", height = 10, width = 20)
for (i in Pheno_list) {
  p_raw <- ggplot(Malay_Phenotype_all_raw, aes(x = .data[[i]])) +
    geom_density(alpha = 0.5, fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 8, face = "bold"))
  
  p_outrem_trans <- ggplot(Malay_Phenotype_all, aes(x = .data[[i]])) +
    geom_density(alpha = 0.5, fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 12, face = "bold"))
  
  combined_plot <- plot_grid(p_raw, p_outrem_trans, labels = c("Raw Data", "Outlier Removed & INT"), ncol = 2)
  
  print(combined_plot)
}
dev.off()


## Split data into multiple files for regenie [each file with FID, IID and 36 Phenotypes]. 
#This will help run Regenie on each file for ease of computation.

IDs <- Malay_Phenotype_all[,c(1,2)]
chunks <- ceiling(1073/36) ## 1073 phenotypes & 36 phenotypes in each file

for (i in 1:chunks) {
  # Calculate the columns for the current chunk
  start_col <- 3 + (i - 1) * 36
  end_col <- min(2 + i * 36, 1075) ## 1075 is the total number of columns
  
  # Select columns for the current chunk
  current_chunk <- as.data.frame(cbind(IDs, Malay_Phenotype_all[, c(start_col:end_col)]))
  
  # Save the current chunk to a file (adjust the filename as needed)
  write.table(current_chunk, file = paste0("Phenotypes/Malay_Phenotype_File_", i, ".txt"), sep = " ", row.names = FALSE,quote = F)
}


####Covariate File for Malay

Malay_Covar <- Covar_file %>% filter(IID %in% Malay_Phenotype_all$IID)
Malay_PCs <- read.table("HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06.eigenvec",header=T)
Malay_Covar <- merge.data.frame(Malay_Covar,Malay_PCs,by="IID")
Malay_Covar <- Malay_Covar[,c(5,1:4,6:24)] ## FID, IID, Age, Sex, Batch, PC1-20
write.table(Malay_Covar,file="Malay_metab_covar.txt",row.names = F,quote = F,sep = " ")

