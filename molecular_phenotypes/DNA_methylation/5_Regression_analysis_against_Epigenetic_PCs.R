### This script uses PRS as an example for running linear regression to calculating its association with the different epigenetic PCs 
### This script is divided into 2 parts.
## Code in this section
##1. Calculate Epigenetic PCs base on the 16444 differentially regulated CpG
##2. Perform regression analysis using PRS as an example

library(dplyr)
library(cowplot)
library(magick)

##1. Calculate Epigenetic PCs base on the 16444 differentially regulated CpG

#load methylation data
load("../beta_QN_rmGSwap_rmDup_HELIOS.RData")

#read file containing the list of 16444 CpG that are differentially expressed between the threee different ethnic group
temp_Markernames_file = "Meth_diff_ethnic_stringent.txt"
Markernames <- read.table(temp_Markernames_file, header = F)
colnames(Markernames)[1] <- "CpG"

cpg_16K = as.character(Markernames$CpG)
beta_16Kcg <- betaRmDup[rownames(betaRmDup) %in% c(cpg_16K),]

#Calculate PCA
df_betaRmDup_pca <- prcomp(t(na.omit(beta_16Kcg)))

#save the PCA results
save(df_betaRmDup_pca, file="df_betaRmDup_pca.RData")

#calculate the eigenvalues base on the top 100 PCs
eigenvals <- df_betaRmDup_pca$sdev^2
eigenvals_top100_por <- eigenvals/sum(eigenvals[c(1:100)])*100

data1 <- data.frame(eigenvals_top100_por)

row.names(data1)[1:5] <- c("PC1", "PC2", "PC3", "PC4", "PC5")
subset_data <- data.frame(
  x = rownames(data1)[1:5],  # Assuming you want row names as x-axis
  y = data1[1:5, 1]           # Assuming you want the first column for the y-axis
)

##2. Perform regression analysis using PRS as an example

sentrixid <- read.csv("../Sentrix_to_HELIOS_NPM_ID_demographic_new.csv")
sentrixid$Gender <- as.factor(sentrixid$Gender)
sentrixid$Ethnicity <- as.factor(sentrixid$Ethnicity)
row.names(sentrixid) <- sentrixid$SentrixID
#str(sentrixid)

#load control probs
load("ctrlprobes.RData")

#load housemen
housemen <- read.table("houseman_constrainedCoefs.txt")
sentrixid_housemen <- merge(sentrixid,housemen, by="row.names")
sentrixid_housemen <- sentrixid_housemen[,-1]
row.names(sentrixid_housemen) <- sentrixid_housemen$SentrixID
sentrixid_housemen_ctrl <-  merge(sentrixid_housemen,ctrlprobes.scores, by = "row.names")
sentrixid_housemen_ctrl <- sentrixid_housemen_ctrl[,-1]

#Read and merge profile PRS

#1
PRS <- read.csv("sentrixid_PRSv2.csv", header = T)
PRSv2 <- PRS[,7:ncol(PRS)]
PRS_temp <- data.frame(PRSv2)
PRSv3 <- cbind(PRS$FREG1_Barcode,PRS_temp)
colnames(PRSv3)[1] <- "HELIOS_ID"

#Merge PRS with sentrixid
sentrixid_PRS <- merge(sentrixid_housemen_ctrl,PRSv3, by = "HELIOS_ID")
sentrixid_PRS[sentrixid_PRS == -Inf] <- NA

sentrixid_PRS_PCA_noothers <- sentrixid_PRS %>%
  subset(Ethnicity != 4)

#start loop to loop through the 5 epigenetic PCs
for (z in 8:12){
  
  #define formula
  lfla_lm_PRS=as.formula("scale(sentrixid_PRS_PCA_noothers[,i]) ~ scale(sentrixid_PRS_PCA_noothers[,z]) + sentrixid_PRS_PCA_noothers$Age + sentrixid_PRS_PCA_noothers$Gender + sentrixid_PRS_PCA_noothers$Ethnicity + sentrixid_PRS_PCA_noothers$CD8T + sentrixid_PRS_PCA_noothers$CD4T + sentrixid_PRS_PCA_noothers$NK + sentrixid_PRS_PCA_noothers$Bcell + sentrixid_PRS_PCA_noothers$Mono + sentrixid_PRS_PCA_noothers$Gran + sentrixid_PRS_PCA_noothers$PC1_cp + sentrixid_PRS_PCA_noothers$PC2_cp + sentrixid_PRS_PCA_noothers$PC3_cp + sentrixid_PRS_PCA_noothers$PC4_cp + sentrixid_PRS_PCA_noothers$PC5_cp + sentrixid_PRS_PCA_noothers$PC6_cp + sentrixid_PRS_PCA_noothers$PC7_cp + sentrixid_PRS_PCA_noothers$PC8_cp + sentrixid_PRS_PCA_noothers$PC9_cp + sentrixid_PRS_PCA_noothers$PC10_cp + sentrixid_PRS_PCA_noothers$PC11_cp + sentrixid_PRS_PCA_noothers$PC12_cp + sentrixid_PRS_PCA_noothers$PC13_cp + sentrixid_PRS_PCA_noothers$PC14_cp + sentrixid_PRS_PCA_noothers$PC15_cp + sentrixid_PRS_PCA_noothers$PC16_cp + sentrixid_PRS_PCA_noothers$PC17_cp + sentrixid_PRS_PCA_noothers$PC18_cp + sentrixid_PRS_PCA_noothers$PC19_cp + sentrixid_PRS_PCA_noothers$PC20_cp + sentrixid_PRS_PCA_noothers$PC21_cp + sentrixid_PRS_PCA_noothers$PC22_cp + sentrixid_PRS_PCA_noothers$PC23_cp + sentrixid_PRS_PCA_noothers$PC24_cp + sentrixid_PRS_PCA_noothers$PC25_cp + sentrixid_PRS_PCA_noothers$PC26_cp + sentrixid_PRS_PCA_noothers$PC27_cp + sentrixid_PRS_PCA_noothers$PC28_cp + sentrixid_PRS_PCA_noothers$PC29_cp + sentrixid_PRS_PCA_noothers$PC30_cp")
  
  #define dimension of martix
  res_lm_PRS=matrix(ncol=4,nrow=ncol(sentrixid_PRS_PCA_noothers[228:ncol(sentrixid_PRS_PCA_noothers)]))
  colnames(res_lm_PRS) <- c('Estimate','Std.Error','z_value','P')
  
  timestamp()

  loop through the different PRS
  for(i in 228:ncol(sentrixid_PRS_PCA_noothers)) {
    tryCatch({fit = summary(lm(lfla_lm_PRS, na.action=na.omit))}, error = function(error) {return(NA)})
    if(!exists("fit")){
      res_lm_PRS[i-227,] = rep(NA,4)
    }else{
      res_lm_PRS[i-227,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
      rm(fit)
    }
  }
  
  timestamp()
  res_lm_PRSv2 <- cbind(colnames(sentrixid_PRS_PCA_noothers)[228:ncol(sentrixid_PRS_PCA_noothers)],res_lm_PRS)
  colnames(res_lm_PRSv2)[1] <- paste0("PRS")
  
  write.csv(res_lm_PRSv2, file=paste0("../PRS/PCA/",colnames(sentrixid_PRS_PCA_noothers)[z],"_results_linear_regression_HELIOS_CpG_PRS.csv"), row.names = F)
  
  rm(res_lm_PRSv2)
  
}

#end of regresison analysis
#______________

