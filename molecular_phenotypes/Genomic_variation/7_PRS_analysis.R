#### Set Path to directory and load packages

setwd("~/path/to/Genomic_Data")

#### Load required packages

library(ggplot2)
library(dplyr)
library(metafor)


##### Code Shown for 1 disorder, same repeated for others

#### PRS Distribution to determine ethnic differences

Chinese_T2D <- read.table("HELIOS_Chinese_T2D.profile",header = T)
Chinese_T2D$Ethnicity <- "Chinese"
Chinese_T2D <- Chinese_T2D[,c(2,7,6)]

Indian_T2D <- read.table("HELIOS_Indian_T2D.profile",header = T)
Indian_T2D$Ethnicity <- "Indian"
Indian_T2D <- Indian_T2D[,c(2,7,6)]

Malay_T2D <- read.table("HELIOS_Malay_T2D.profile",header = T)
Malay_T2D$Ethnicity <- "Malay"
Malay_T2D <- Malay_T2D[,c(2,7,6)]

T2D_all <- as.data.frame(rbind(Chinese_T2D,Indian_T2D,Malay_T2D))

T2D_all$SCORE <- (T2D_all$SCORE - mean(T2D_all$SCORE))/(sd(T2D_all$SCORE))




Phenotype_info <- read.table("HELIOS_T2D_status.txt",header = F) ## refer to descriptive_health_outcomes_script_for_information
names(Phenotype_info) <- c("IID","T2D_Status")

Covar_file <- read.table("HELIOS_Covar_file.txt",header = T) ### Age and Gender Information
names(Covar_file) <- c("IID","Age","Gender")

All_Info <- merge.data.frame(Phenotype_info,Covar_file,by="IID")


#### ANOVA test to determine significant difference between ethnic groups

T2D_all <- merge.data.frame(T2D_all,All_Info,by="IID")
T2D_anova <- aov(SCORE~Ethnicity+Age+Gender,T2D_all)
p_T2D <- summary(T2D_anova)[[1]][[5]][[1]]




###### Odds ratio to determine the genetic risk of T2D and the predictive power

Chinese_Final <- merge.data.frame(Chinese_T2D,All_Info,by="IID")
T2D_Chinese<- glm(T2D~1+SCORE+Age+Gender,data=Chinese_Final,family = "binomial")
T2D_Chinese_OR <- exp(coef(T2D_Chinese)[[2]])
T2D_Chinese_OR_lower <- exp(confint(T2D_Chinese,level = 0.95)[[2]])
T2D_Chinese_OR_upper <- exp(confint(T2D_Chinese,level = 0.95)[[6]])
T2D_Chinese_P <- summary(T2D_Chinese)$coefficients[, "Pr(>|z|)"][[2]]

Indian_Final <- merge.data.frame(Indian_T2D,All_Info,by="IID")
T2D_Indian<- glm(T2D~1+SCORE+Age+Gender,data=Indian_Final,family = "binomial")
T2D_Indian_OR <- exp(coef(T2D_Indian)[[2]])
T2D_Indian_OR_lower <- exp(confint(T2D_Indian,level = 0.95)[[2]])
T2D_Indian_OR_upper <- exp(confint(T2D_Indian,level = 0.95)[[6]])
T2D_Indian_P <- summary(T2D_Indian)$coefficients[, "Pr(>|z|)"][[2]]

Malay_Final <- merge.data.frame(Malay_T2D,All_Info,by="IID")
T2D_Malay<- glm(T2D~1+SCORE+Age+Gender,data=Malay_Final,family = "binomial")
T2D_Malay_OR <- exp(coef(T2D_Malay)[[2]])
T2D_Malay_OR_lower <- exp(confint(T2D_Malay,level = 0.95)[[2]])
T2D_Malay_OR_upper <- exp(confint(T2D_Malay,level = 0.95)[[6]])
T2D_Malay_P <- summary(T2D_Malay)$coefficients[, "Pr(>|z|)"][[2]]



#### Trans-ancestry meta analysis and Heterogeneity


T2D_beta <- c(coefficients(T2D_Chinese)[[2]],
              coefficients(T2D_Indian)[[2]],
              coefficients(T2D_Malay)[[2]])

T2D_se <- c(summary(T2D_Chinese)$coefficients[,2][[2]],
            summary(T2D_Indian)$coefficients[,2][[2]],
            summary(T2D_Malay)$coefficients[,2][[2]])
T2D_data <- data.frame(beta = T2D_beta, se = T2D_se)

T2D_meta <- rma(yi = T2D_data$bet, sei = T2D_data$se)
summary(T2D_meta)
T2D_trans_OR <- as.vector(exp(T2D_meta$b))
T2D_trans_P <- as.vector(T2D_meta$pval)
T2D_trans_lower <- as.vector(exp(T2D_meta$ci.lb))
T2D_trans_upper <- as.vector(exp(T2D_meta$ci.ub))
T2D_I2 <- T2D_meta$I2
T2D_hetp<- T2D_meta$QEp



###########
