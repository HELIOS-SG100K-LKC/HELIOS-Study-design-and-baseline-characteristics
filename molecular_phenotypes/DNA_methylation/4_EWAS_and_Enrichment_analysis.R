### This script uses PRS as an example for running EWAS and calculate the percentage enrichment of ethnic different CgG that are associated with the different PRS.
### This script is divided into 2 parts.
## Code in this section
#1. Run EWAS code
#2. Extract ethnic different CpG and Perforem hypergeometric test (phyper) and calculate enrichment


#1. Run EWAS code
library(dplyr)

#oad methylation data
load("beta_QN_rmGSwap_rmDup_HELIOS.RData")

#load control probs
load("ctrlprobes.RData")

#load housemen
housemen <- read.table("houseman_constrainedCoefs.txt")

#read sentrix id, helios id and merge helios id with contorl probes and housemen data
sentrixid <- read.csv("Sentrix_to_HELIOS_NPM_ID_demographic.csv", header = T)
colnames(sentrixid)[2] <- "FREG1_Barcode"
sentrixid$Gender <- as.factor(sentrixid$Gender)
sentrixid$Ethnicity <- as.factor(sentrixid$Ethnicity)
row.names(sentrixid) <- sentrixid$SentrixID

sentrixid_housemen <- merge(sentrixid,housemen, by="row.names")
sentrixid_housemen <- sentrixid_housemen[,-1]
row.names(sentrixid_housemen) <- sentrixid_housemen$SentrixID
sentrixid_housemen_ctrl <-  merge(sentrixid_housemen,ctrlprobes.scores, by = "row.names")
sentrixid_housemen_ctrl <- sentrixid_housemen_ctrl[,-1]

#Prepare phenotype files
#read PRS
PRS <- read.csv("../PRS/sentrixid_PRSv2.csv", header = T)
PRSv2 <- PRS[,7:ncol(PRS)]
PRS_temp <- data.frame(PRSv2)
PRSv3 <- cbind(PRS$FREG1_Barcode,PRS_temp)
colnames(PRSv3)[1] <- "FREG1_Barcode"
PRSv3[PRSv3 == -Inf] <- NA

#merge sentrixid file with phenotype files
sentrixid_PRSv3_housemen_ctrl <- merge(sentrixid_housemen_ctrl,PRSv3, by="FREG1_Barcode")

#filter out participants where ethnicity is not Chinese, Malat or Indian
sentrixid_PRSv3_housemen_ctrl_noothers <- sentrixid_PRSv3_housemen_ctrl %>% filter(Ethnicity != 4)

#sort by sentrixid
sentrixid_PRSv3_housemen_ctrl_noothers_sorted <- sentrixid_PRSv3_housemen_ctrl_noothers[order(sentrixid_PRSv3_housemen_ctrl_noothers$SentrixID),]

#extract from betaRmDup particiapnts that we and analysing and sort column name which is the sentrixid
samples <- sentrixid_PRSv3_housemen_ctrl_noothers$SentrixID
betaRmDup_PRSv3 <- betaRmDup[,colnames(betaRmDup) %in% samples]
betaRmDup_PRSv3_sorted <- betaRmDup_PRSv3[,order(as.character(colnames(betaRmDup_PRSv3)))]

#perform linear regression

#start loop
for (z in c(217:221)){

  print(z)
  
  #define model
  lfla_PRS=as.formula("sentrixid_PRSv3_housemen_ctrl_noothers_sorted[,z] ~ scale(betaRmDup_PRSv3_sorted[i,]) + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$Age + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$Gender + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$Ethnicity + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$CD8T + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$CD4T + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$NK + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$Bcell + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$Mono + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$Gran + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC1_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC2_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC3_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC4_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC5_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC6_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC7_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC8_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC9_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC10_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC11_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC12_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC13_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC14_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC15_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC16_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC17_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC18_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC19_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC20_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC21_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC22_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC23_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC24_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC25_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC26_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC27_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC28_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC29_cp + sentrixid_PRSv3_housemen_ctrl_noothers_sorted$PC30_cp")
  
  nvar = nrow(betaRmDup_PRSv3_sorted)
  HELIOS_EWAS_PRS=matrix(ncol=4, nrow=nvar)

  timestamp()

  colnames(HELIOS_EWAS_PRS) =c('Estimate','Std.Error', 'z_value','P')  

  
for(i in 1:nvar) {
  tryCatch({fit = summary(lm(lfla_PRS))}, error = function(error) {return(NA)})
  if(!exists("fit")){
    HELIOS_EWAS_PRS[i,] = rep(NA,4)
  }else{
    HELIOS_EWAS_PRS[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
    rm(fit)
  }
}

timestamp()
#HELIOS_EWAS_PRS$CpG <- rownames(betaRmDup_PRSv3_sorted)
rownames(HELIOS_EWAS_PRS) <- rownames(betaRmDup_PRSv3_sorted)

write.csv(HELIOS_EWAS_PRS, file=paste0("../PRS/",colnames(sentrixid_PRSv3_housemen_ctrl_noothers_sorted)[z],"_HELIOS_EWAS_PRS.csv"), row.names = T)
save(HELIOS_EWAS_PRS, file=paste0("../PRS/",colnames(sentrixid_PRSv3_housemen_ctrl_noothers_sorted)[z],"_HELIOS_EWAS_PRS.RData"))

}

#2. Extract ethnic different CpG and Perforem hypergeometric test (phyper) and calculate enrichment

##Prepare a file that contain the list of all the traits and group that we want to calcuate the enrichment for
#read variables
variable_list <- read.csv("EWAS_Grp_variable_Methods.csv", header = T)

#read the file that contains the list of 16444 CpG that are differentially regulated in the 3 different ethnic group
eff <- read.csv("Meth_diff_ethnic_stringent.txt", header = F)
colnames(eff)[1] <- "CpG"

#Begin extraction
Enrichment_16K_summary <- data.frame()

for (i in 1:nrow(variable_list)){

  print(i)
  
  timestamp()

grp <- variable_list[i,1]
variable <- variable_list[i,2]

#temp_Markernames_file = "../../Sentinel_and_PCs_CpG_regression_phewas_PRS/Meth_diff_ethnic_stringent.txt"
#Markernames <- read.table(temp_Markernames_file, header = T)
#colnames(Markernames)[1] <- "CpG"
cpg_16K = as.character(eff$CpG)

#read EWAS results
res_EWAS <- read.csv(paste0("../",grp,"/",variable,"_HELIOS_EWAS_",grp,".csv"))

#Extract the results from the list of 16444 CpG of interest  
res_16K <- res_EWAS[res_EWAS$X %in% c(cpg_16K),]
colnames(res_16K)[1] <- "CpG"

#write the results of the 16444 Cpg into a different file for future reference
write.csv(res_16K,file=paste0("../../16K_Cpg/",variable,"_HELIOS_16K_Cpgs_",grp,".csv"), row.names = F)


#Perform hypergeometric testing using phyper  
P = nrow(res_EWAS)
#X = 5000
N = nrow(res_16K)
#x = 300

x = sum(res_16K$P < 0.05)
X = sum(res_EWAS$P < 0.05)

sig_EWAS = X/P
enrich_16K = x/N

Phyper_16K = as.numeric(phyper(x-1,X,P-X,N,lower.tail=F))
fold_enrich = as.numeric(enrich_16K/sig_EWAS)
Enrichment_16K = cbind(grp,variable,x,N,enrich_16K,X,P,sig_EWAS,fold_enrich,Phyper_16K)

Enrichment_16K_summary <- rbind(Enrichment_16K_summary,Enrichment_16K)

rm(Enrichment_16K)
rm(res_16K)
rm(res_EWAS)

timestamp()

}

#write the results out into a file
write.csv(Enrichment_16K_summary,file="../../16K_Cpg/Enrichment_16K_summary.csv")


