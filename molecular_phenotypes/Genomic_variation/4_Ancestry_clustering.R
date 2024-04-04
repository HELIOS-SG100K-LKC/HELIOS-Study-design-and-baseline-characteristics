#### Set Path to directory and load packages

setwd("~/path/to/Genomic_Data")

library(ggplot2)
library(dplyr)
library(cluster)
library(cowplot)
library(tidyr)
library(ggstar)
library(ggbreak)

#### load and prepare files
PCs <- read.table("HELIOS_All.eigenvec",header=T)

Ancestry_data <- read.delim("HELIOS_reported_Ancestry_data.txt",header = T)
names(Ancestry_data) <- c("FID","IID","Reported_Ancestry")
Ancestry_data <- Ancestry_data %>% filter(IID%in% PCs$IID)

Data_all <- merge.data.frame(PCs,Ancestry_data,by=c("FID","IID"))


##### Inital K- means clustering [To select homogenous individuals from 3 main ancestry groups]
selected_vars <- Data_all[,c(12:16)]

kmeans_model_1 <- kmeans(selected_vars, centers = 3, nstart = 10)
Data_all$cluster <- kmeans_model_1$cluster
print(kmeans_model_1$centers)

## Clusters and reported ancesry matching [1: Chinese; 2: Malay; 3: Indian ]


centers <- kmeans_model_1$centers
centers <- as.data.frame(centers)
rownames(centers) <- c("Chinese","Malay","Indian")

## Select Individuals within 2 SD of Chinese PC 1-5 Centers

sd_pc1 <- Data_all %>% filter(cluster2==1) %>% select(PC1) %>% summarise(sd=sd(PC1)) %>% pull(sd)
sd_pc2 <- Data_all %>% filter(cluster2==1) %>% select(PC2) %>% summarise(sd=sd(PC2)) %>% pull(sd)
sd_pc3 <- Data_all %>% filter(cluster2==1) %>% select(PC3) %>% summarise(sd=sd(PC3)) %>% pull(sd)
sd_pc4 <- Data_all %>% filter(cluster2==1) %>% select(PC4) %>% summarise(sd=sd(PC4)) %>% pull(sd)
sd_pc5 <- Data_all %>% filter(cluster2==1) %>% select(PC5) %>% summarise(sd=sd(PC5)) %>% pull(sd)

pc1_upper <- centers["Chinese","PC1"] + 2*sd_pc1
pc1_lower <- centers["Chinese","PC1"] - 2*sd_pc1
pc2_upper <- centers["Chinese","PC2"] + 2*sd_pc2
pc2_lower <- centers["Chinese","PC2"] - 2*sd_pc2
pc3_upper <- centers["Chinese","PC3"] + 2*sd_pc3
pc3_lower <- centers["Chinese","PC3"] - 2*sd_pc3
pc4_upper <- centers["Chinese","PC4"] + 2*sd_pc4
pc4_lower <- centers["Chinese","PC4"] - 2*sd_pc4
pc5_upper <- centers["Chinese","PC5"] + 2*sd_pc5
pc5_lower <- centers["Chinese","PC5"] - 2*sd_pc5

Data_all$anc_group_PC1_5 <- ifelse(Data_all$PC1<=pc1_upper & Data_all$PC1>=pc1_lower & 
                                     Data_all$PC2<=pc2_upper & Data_all$PC2>=pc2_lower & 
                                     Data_all$PC3<=pc3_upper & Data_all$PC3>=pc3_lower & 
                                     Data_all$PC4<=pc4_upper & Data_all$PC4>=pc4_lower & 
                                     Data_all$PC5<=pc5_upper & Data_all$PC5>=pc5_lower, "Chinese","Others")


## Select Individuals within 2 SD of Malay PC 1-5 Centers

sd_pc1 <- Data_all %>% filter(cluster2==2) %>% select(PC1) %>% summarise(sd=sd(PC1)) %>% pull(sd)
sd_pc2 <- Data_all %>% filter(cluster2==2) %>% select(PC2) %>% summarise(sd=sd(PC2)) %>% pull(sd)
sd_pc3 <- Data_all %>% filter(cluster2==2) %>% select(PC3) %>% summarise(sd=sd(PC3)) %>% pull(sd)
sd_pc4 <- Data_all %>% filter(cluster2==2) %>% select(PC4) %>% summarise(sd=sd(PC4)) %>% pull(sd)
sd_pc5 <- Data_all %>% filter(cluster2==2) %>% select(PC5) %>% summarise(sd=sd(PC5)) %>% pull(sd)

pc1_upper <- centers["Malay","PC1"] + 2*sd_pc1
pc1_lower <- centers["Malay","PC1"] - 2*sd_pc1
pc2_upper <- centers["Malay","PC2"] + 2*sd_pc2
pc2_lower <- centers["Malay","PC2"] - 2*sd_pc2
pc3_upper <- centers["Malay","PC3"] + 2*sd_pc3
pc3_lower <- centers["Malay","PC3"] - 2*sd_pc3
pc4_upper <- centers["Malay","PC4"] + 2*sd_pc4
pc4_lower <- centers["Malay","PC4"] - 2*sd_pc4
pc5_upper <- centers["Malay","PC5"] + 2*sd_pc5
pc5_lower <- centers["Malay","PC5"] - 2*sd_pc5


Data_all$anc_group_PC1_5 <- ifelse(Data_all$PC1<=pc1_upper & Data_all$PC1>=pc1_lower & 
                                     Data_all$PC2<=pc2_upper & Data_all$PC2>=pc2_lower & 
                                     Data_all$PC3<=pc3_upper & Data_all$PC3>=pc3_lower & 
                                     Data_all$PC4<=pc4_upper & Data_all$PC4>=pc4_lower & 
                                     Data_all$PC5<=pc5_upper & Data_all$PC5>=pc5_lower, "Malay",Data_all$anc_group_PC1_5)



## Select Individuals within 2 SD of Indian PC 1-5 Centers
sd_pc1 <- Data_all %>% filter(cluster2==3) %>% select(PC1) %>% summarise(sd=sd(PC1)) %>% pull(sd)
sd_pc2 <- Data_all %>% filter(cluster2==3) %>% select(PC2) %>% summarise(sd=sd(PC2)) %>% pull(sd)
sd_pc3 <- Data_all %>% filter(cluster2==3) %>% select(PC3) %>% summarise(sd=sd(PC3)) %>% pull(sd)
sd_pc4 <- Data_all %>% filter(cluster2==3) %>% select(PC4) %>% summarise(sd=sd(PC4)) %>% pull(sd)
sd_pc5 <- Data_all %>% filter(cluster2==3) %>% select(PC5) %>% summarise(sd=sd(PC5)) %>% pull(sd)

pc1_upper <- centers["Indian","PC1"] + 2*sd_pc1
pc1_lower <- centers["Indian","PC1"] - 2*sd_pc1
pc2_upper <- centers["Indian","PC2"] + 2*sd_pc2
pc2_lower <- centers["Indian","PC2"] - 2*sd_pc2
pc3_upper <- centers["Indian","PC3"] + 2*sd_pc3
pc3_lower <- centers["Indian","PC3"] - 2*sd_pc3
pc4_upper <- centers["Indian","PC4"] + 2*sd_pc4
pc4_lower <- centers["Indian","PC4"] - 2*sd_pc4
pc5_upper <- centers["Indian","PC5"] + 2*sd_pc5
pc5_lower <- centers["Indian","PC5"] - 2*sd_pc5


Data_all$anc_group_PC1_5 <- ifelse(Data_all$PC1<=pc1_upper & Data_all$PC1>=pc1_lower & 
                                     Data_all$PC2<=pc2_upper & Data_all$PC2>=pc2_lower & 
                                     Data_all$PC3<=pc3_upper & Data_all$PC3>=pc3_lower & 
                                     Data_all$PC4<=pc4_upper & Data_all$PC4>=pc4_lower & 
                                     Data_all$PC5<=pc5_upper & Data_all$PC5>=pc5_lower, "Indian",Data_all$anc_group_PC1_5)



## Get Count of Individuals within each group
Data_all %>% group_by(anc_group_PC1_5) %>% summarise(count=n())

## Counts [Chinese: 5414; Indian: 1141; Malay: 981; Others:2401]


#### Write cluster Info to perform supervised Admixture analysis using SCOPE

Ancestry_Info_clust3 <- Data_all[,c(1,2,24)]
names(Ancestry_Info_clust3) <- c("FID","IID","Cluster")
write.table(Ancestry_Info_clust3,file="Homogenous_ancestry_info.txt",row.names = F,quote = F,sep = " ")

############# RUN SCOPE using run_SCOPE.sh; use output of SCOPE for remaining steps ###################3


#### Final Clustering Using Admxiture Proportions

Admix_prop <- read.table("HELIOS_All_Supervised_PC1_5.Qhat.txt",header = F)
Admix_prop <- as.data.frame(t(Admix_prop))
Admix_prop$IID <- PCs$IID
names(Admix_prop) <- c("Pop_1","Pop_2","Pop_3","IID")

Data_all <- merge.data.frame(Data_all,Admix_prop,by="IID")



#### K- means clustering using Admixture proportions (K=3, no defined centers)
selected_vars <- Data_all[,c(25:27)]


kmeans_model_2 <- kmeans(standardized_vars, centers = 3, nstart = 10)


#### get centers
cet_matrix <- kmeans_model_2$centers
cet_matrix <- as.data.frame(cet_matrix)
cet_matrix$ancestry <- c("C","I","M")

#### Define Centers for Admixed Ancestry

value_CI <- colMeans(cet_matrix[c(1,2), c(1:3)])
value_CI <- c(value_CI,"CI")
value_CM <- colMeans(cet_matrix[c(1,3), c(1:3)])
value_CM <- c(value_CM,"CM")
value_IM <- colMeans(cet_matrix[c(2,3), c(1:3)])
value_IM <- c(value_IM,"IM")


cet_matrix <- as.data.frame(rbind(cet_matrix,value_CI,value_CM,value_IM))
cet_matrix$Pop_1 <- as.numeric(as.vector(cet_matrix$Pop_1))
cet_matrix$Pop_2 <- as.numeric(as.vector(cet_matrix$Pop_2))
cet_matrix$Pop_3 <- as.numeric(as.vector(cet_matrix$Pop_3))

centroids <- as.matrix(cet_matrix[,c(1:3)])

### Final K - means clustering using the 6 defined centers
kmeans_model_fimnal <- kmeans(standardized_vars, centers = centroids, nstart = 10,iter.max = 10)
Data_all$cluster_final <- kmeans_model_final$cluster

#kmeans_model2$centers

## get count for each group
Data_all %>% group_by(cluster_k6) %>% summarise(count=n())

## Count = [Chinese: 7079; Indian: 1507; Malay: 852; Chinese-Indian: 115; Chinese-Malay: 203; Indian-Malay: 181]


#### Extract and Save Final clustering Information

Final_Cluster_Info <- Data_all[,c(1,2,28)]

Pred_anc <- as.data.frame(cbind(c(1:6),c("Chinese","Indian","Malay","Chinese-Indian","Chinese-Malay","Indian-Malay")))
names(Pred_anc) <- c("cluster_k6","Predicted_Ancestry")

Final_Cluster_Info <- merge.data.frame(Final_Cluster_Info,Pred_anc,by="cluster_k6")
write.table(Final_Cluster_Info,file = "Final_Clustering_Information.txt",row.names = F,quote = F,sep=" ")



#####
