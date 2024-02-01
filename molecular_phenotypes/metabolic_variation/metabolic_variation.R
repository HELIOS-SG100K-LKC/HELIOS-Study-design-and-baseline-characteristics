
## This script has 6 parts (Step 1-6): 
    # Step 1: GLM for metabolite-ethnicity associations
    # Step 2: Validate significant results in test set 
    # Step 3: PCA of 153 metabolites 
    # Step 4: Correlation matrix for metabolites 
    # Step 5: Generate Circos plot 
    # Step 6: Metabolic variance explained

## NOTE: 
  # Results from Steps 1-2 are made available and is used in consequent analyses in Steps 3-6.
  # Results from Step 4 are made available and used in Step 5. 
  # Aggregate results of association of metabolites with HT, Obs, T2D, and CVD are also made available and used in Step 5. 

## This script was run on R4.2.1



#####################
## load libraries ##
#####################

library(data.table)
library(tidyr)
library(dplyr)
library(ggfortify)
library(glue)
library(stringr)
library(tibble)
library(forcats)
library(Hmisc)
library(circlize)
library(ComplexHeatmap)
library(ppcor)



########################################################
## Step 1: GLM for metabolite-ethnicity associations ##
########################################################

## remove environment variables
rm(list = ls())


## read phenotype file
df_pheno <- fread("<insert phenotype filename>", 
                  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


## format phenotype data for analysis
# remove Ethnicity 'others'
df_pheno %>%
  filter(Ethnicity != "O") -> df_pheno

# subset and factorize phenotypes 
df_pheno$Ethnicity = as.factor(df_pheno$Ethnicity)
df_pheno$Sex = as.factor(df_pheno$Sex)
df_pheno$SSID = as.factor(df_pheno$SSID) #SSID is batch number

# drop na 
df_pheno %>% 
  drop_na(Age, Sex, Ethnicity) -> df_pheno


## divide phenotype data into test and discovery cohorts
set.seed(1234)
train_size <- floor(0.7*nrow(df_pheno))
trainrows <- sample(1:nrow(df_pheno), size = train_size, replace = FALSE)

df_pheno_train <- df_pheno[trainrows, ]
df_pheno_train %>% arrange(PARENT_SAMPLE_NAME) -> df_pheno_train

df_pheno_test <- df_pheno[-trainrows, ]
df_pheno_test %>% arrange(PARENT_SAMPLE_NAME) -> df_pheno_test


## read metabolite data
df_metab <- fread("<insert metabolite filename>", 
            header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


## create train and test data for metabolite and sort samples in same order as phenotype data
df_metab %>% 
  filter(PARENT_SAMPLE_NAME %in% df_pheno_train$PARENT_SAMPLE_NAME) %>% 
  arrange(PARENT_SAMPLE_NAME) -> df_metab_train

df_metab %>% 
  filter(PARENT_SAMPLE_NAME %in% df_pheno_test$PARENT_SAMPLE_NAME) %>% 
  arrange(PARENT_SAMPLE_NAME) -> df_metab_test

table(df_metab_train$PARENT_SAMPLE_NAME == df_pheno_train$PARENT_SAMPLE_NAME) #all true only if in same order
table(df_metab_test$PARENT_SAMPLE_NAME == df_pheno_test$PARENT_SAMPLE_NAME) #all true only if in same order


## format metabolite train set for analysis
# retain all columns except sample names
df_metab_train = df_metab_train[, -1] 

# remove metabolites with column variance = 0 
df_metab_train %>% 
  dplyr::select_if(~var(., na.rm = TRUE) != 0) -> df_metab_train

# save metabolite names
CHEM_ID <- colnames(df_metab_train)

# convert to matrix
colnames(df_metab_train) = NULL
df_metab_train <- as.matrix(df_metab_train)


## run linear regression
# NOTE: Chinese automatically set as ethnicity reference 
lm.formula <- as.formula('df_metab_train[, i] ~ df_pheno_train$Ethnicity + df_pheno_train$Age + df_pheno_train$Sex + df_pheno_train$SSID')

res_I = data.frame(matrix(nrow = ncol(df_metab_train), ncol = 4))
colnames(res_I) =c('Est','SE', 't_value','P')

res_M = data.frame(matrix(nrow = ncol(df_metab_train), ncol = 4))
colnames(res_M) =c('Est','SE', 't_value','P')

system.time(
  for(i in 1:ncol(df_metab_train)) {
    tryCatch({reg.out <- summary(glm(lm.formula, family = gaussian(), na.action = na.omit))}, error = function(error) {return(NA)})
    
    if(!exists("reg.out")) {
      res_I[i,] = res_M[i,] = rep(NA,4)
    } else {
      res_I[i,] = tryCatch({reg.out$coefficients[2,]},error = function(error) {return(rep(NA,4))})
      res_M[i,] = tryCatch({reg.out$coefficients[3,]},error = function(error) {return(rep(NA,4))})
      
      rm(reg.out)
    }
    if(i%%100 == 0) {
      print(i)
    }
  }
)

res_I_final <- cbind(CHEM_ID, res_I, stringsAsFactors = FALSE)
res_M_final <- cbind(CHEM_ID, res_M, stringsAsFactors = FALSE)

head(res_I_final)
head(res_M_final)

res_I_final %>% 
  dplyr::rename_with(~paste0(., "_I"), 2:5) -> res_I_final

res_M_final %>% 
  dplyr::rename_with(~paste0(., "_M"), 2:5) -> res_M_final


## For Indian vs Malay comparison set ethnicity reference as Indian
df_pheno_train <- within(df_pheno_train, Ethnicity <- relevel(Ethnicity, ref = "I"))

lm.formula <- as.formula('df_metab_train[, i] ~ df_pheno_train$Ethnicity + df_pheno_train$Age + df_pheno_train$Sex + df_pheno_train$SSID')

res_IM = data.frame(matrix(nrow = ncol(df_metab_train), ncol = 4))
colnames(res_IM) =c('Est','SE', 't_value','P')

system.time(
  for(i in 1:ncol(df_metab_train)) {
    tryCatch({reg.out <- summary(glm(lm.formula, family = gaussian(), na.action = na.omit))}, error = function(error) {return(NA)})
    
    if(!exists("reg.out")) {
      res_IM[i,] = rep(NA,4)
    } else {
      res_IM[i,] = tryCatch({reg.out$coefficients[3,]},error = function(error) {return(rep(NA,4))})
      
      rm(reg.out)
    }
    if(i%%100 == 0) {
      print(i)
    }
  }
)

res_IM_final <- cbind(CHEM_ID, res_IM, stringsAsFactors = FALSE)
head(res_IM_final)

res_IM_final %>% 
  dplyr::rename_with(~paste0(., "_IM"), 2:5) -> res_IM_final


## merge all results
res_I_final %>% 
  left_join(res_M_final, by = "CHEM_ID") %>% 
  left_join(res_IM_final, by = "CHEM_ID") -> train_out

## get chem annotations
df_chem <- fread("<insert chemical annotation filename>", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)

df_chem %>%
  mutate_if(is.character, str_trim) -> df_chem
df_chem$CHEM_ID = as.character(df_chem$CHEM_ID)

## retrieve significant results with annotations
train_out %>% 
  left_join(df_chem, by = "CHEM_ID") %>% 
  dplyr::select(SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, colnames(train_out)) %>%
  filter(P_I < 1e-05 & P_M < 1e-05 & P_IM < 1e-05) %>% #1073*3 tests, bonf-p=1e-05
  arrange(SUPER_PATHWAY) -> train_out_sig
#there are 162 significant results



######################################################
## Step 2: Validate significant results in test set ## 
######################################################

## format metabolite test set for analysis
# subset test dataset for significant metabolites
df_metab_test %>% 
  dplyr::select(which(colnames(df_metab_test) %in% train_out_sig$CHEM_ID)) -> df_metab_test 

# remove metabolites with column variance = 0
df_metab_test %>% 
  dplyr::select_if(~var(., na.rm = TRUE) != 0) -> df_metab_test

# save metabolite names
CHEM_ID <- colnames(df_metab_test)

# convert to matrix
colnames(df_metab_test) = NULL
df_metab_test <- as.matrix(df_metab_test)


## run linear regression
# NOTE: Chinese automatically set as ethnicity reference 
lm.formula <- as.formula('df_metab_test[, i] ~ df_pheno_test$Ethnicity + df_pheno_test$Age + df_pheno_test$Sex + df_pheno_test$SSID')

res_I = data.frame(matrix(nrow = ncol(df_metab_test), ncol = 4))
colnames(res_I) =c('Est','SE', 't_value','P')

res_M = data.frame(matrix(nrow = ncol(df_metab_test), ncol = 4))
colnames(res_M) =c('Est','SE', 't_value','P')

system.time(
  for(i in 1:ncol(df_metab_test)) {
    tryCatch({reg.out <- summary(glm(lm.formula, family = gaussian(), na.action = na.omit))}, error = function(error) {return(NA)})
    
    if(!exists("reg.out")) {
      res_I[i,] = res_M[i,] = rep(NA,4)
    } else {
      res_I[i,] = tryCatch({reg.out$coefficients[2,]},error = function(error) {return(rep(NA,4))})
      res_M[i,] = tryCatch({reg.out$coefficients[3,]},error = function(error) {return(rep(NA,4))})
      
      rm(reg.out)
    }
    if(i%%100 == 0) {
      print(i)
    }
  }
)

res_I_final <- cbind(CHEM_ID, res_I, stringsAsFactors = FALSE)
res_M_final <- cbind(CHEM_ID, res_M, stringsAsFactors = FALSE)

head(res_I_final)
head(res_M_final)

res_I_final %>% 
  dplyr::rename_with(~paste0(., "_Itest"), 2:5) -> res_I_final

res_M_final %>% 
  dplyr::rename_with(~paste0(., "_Mtest"), 2:5) -> res_M_final


## For Indian vs Malay comparison set ethnicity reference as Indian
df_pheno_test <- within(df_pheno_test, Ethnicity <- relevel(Ethnicity, ref = "I"))

lm.formula <- as.formula('df_metab_test[, i] ~ df_pheno_test$Ethnicity + df_pheno_test$Age + df_pheno_test$Sex + df_pheno_test$SSID')

res_IM = data.frame(matrix(nrow = ncol(df_metab_test), ncol = 4))
colnames(res_IM) =c('Est','SE', 't_value','P')

system.time(
  for(i in 1:ncol(df_metab_test)) {
    tryCatch({reg.out <- summary(glm(lm.formula, family = gaussian(), na.action = na.omit))}, error = function(error) {return(NA)})
    
    if(!exists("reg.out")) {
      res_IM[i,] = rep(NA,4)
    } else {
      res_IM[i,] = tryCatch({reg.out$coefficients[3,]},error = function(error) {return(rep(NA,4))})
      
      rm(reg.out)
    }
    if(i%%100 == 0) {
      print(i)
    }
  }
)

res_IM_final <- cbind(CHEM_ID, res_IM, stringsAsFactors = FALSE)
head(res_IM_final)

res_IM_final %>% 
  dplyr::rename_with(~paste0(., "_IMtest"), 2:5) -> res_IM_final


## merge all results
res_I_final %>% 
  left_join(res_M_final, by = "CHEM_ID") %>% 
  left_join(res_IM_final, by = "CHEM_ID") -> test_out


## extract metabolites that remain significant and have same direction of effect
test_out %>% 
  filter(P_Itest < 0.05 & P_Mtest < 0.05 & P_IMtest < 0.05) %>% 
  left_join(train_out_sig, by = "CHEM_ID") %>% 
  filter(sign(Est_Itest) == sign(Est_I) & sign(Est_Mtest) == sign(Est_M) & sign(Est_IMtest) == sign(Est_IM)) -> final_out_sig
#153 metabolites are significant in test set


## results available as "ethnicity_metabolites_153_sigall3_testtrain.csv"



#####################################
## Step 3: PCA of 153 metabolites ##
####################################

## subset metabolite data for the 153 significant results
df_metab %>%
  dplyr::select("PARENT_SAMPLE_NAME", which(colnames(df_metab) %in% final_out_sig$CHEM_ID)) -> df_metab_153


## read age-sex matched ids for three ethnicities
matched_ids <- fread("<insert filename of age-sex matched ids>", 
                     header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


## subset matched ids in phenotype and metabolite data
df_pheno %>% 
  filter(`HELIOS Participant ID` %in% matched_ids$`HELIOS Participant ID`) -> df_pheno_matched

df_metab_153 %>% 
  filter(PARENT_SAMPLE_NAME %in% df_pheno_matched$PARENT_SAMPLE_NAME) -> df_metab_153_matched

table(df_metab_153_matched$PARENT_SAMPLE_NAME == df_pheno_matched$PARENT_SAMPLE_NAME) #all true only if in same order


## remove sample names
df_metab_153_matched %>% 
  column_to_rownames("PARENT_SAMPLE_NAME") -> df_forPCA


## remove metabolites with column variance = 0 or with NAs
df_forPCA %>%
  dplyr::select_if(~sum(is.na(.)) == 0) %>% 
  dplyr::select_if(~var(., na.rm = TRUE) != 0) -> df_forPCA


## perform PCA
pca_out <- prcomp(df_forPCA, center = T, scale. = T) 


## generate PCA plot
pca_plot <- as.data.frame(pca_out$x[,1:5])
pca_plot %>% 
  rownames_to_column(var = "PARENT_SAMPLE_NAME") %>%
  left_join(df_pheno_matched, by = "PARENT_SAMPLE_NAME") -> pca_plot

ggplot(pca_plot, aes(x = PC1, y = PC2, col = Ethnicity))+
  geom_point(size = 1.5, alpha = 0.8, aes(shape = Ethnicity))+ 
  scale_color_manual(values = c("#fc8d62", "#8da0cb", "#66c2a5"))+ 
  scale_shape_manual(values = c(19, 15, 17)) +
  stat_ellipse(type = "t", linetype = 1, lwd = 1.2) + 
  theme_classic(18)



#################################################
## Step 4: Correlation matrix for metabolites ##
################################################

## generate correlation matrix for 153 metabolites and format results 
# correlation analysis
rescorr <- Hmisc::rcorr(as.matrix(df_metab_153))

# exttract correlation coefficients
cormat <- rescorr$r
# extract p-values
pmat <- rescorr$P 

# add column and row names
colnames(cormat) <- colnames(df_metab_153)
rownames(cormat) <- colnames(df_metab_153)
colnames(pmat) <- colnames(df_metab_153)
rownames(pmat) <- colnames(df_metab_153)

# convert to long format and merge 
res_r <- reshape2::melt(cormat, value.name = "corr")
res_p <- reshape2::melt(pmat, value.name = "p")
left_join(res_r, res_p, by = c("Var1", "Var2")) -> corr_metab


## subset for significant results
corr_metab %>% 
  drop_na() %>% 
  filter(p < 1e-06) -> corr_metab


## get chem annotations
df_chem <- fread("<insert chemical annotation filename>", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)

df_chem %>%
  mutate_if(is.character, str_trim) -> df_chem
df_chem$CHEM_ID = as.numeric(df_chem$CHEM_ID)


## add chemical annotations for all metabolites in the correlation results
df_chem %>% 
  mutate(Var1 = CHEM_ID) %>% 
  dplyr::select(Var1, CHEMICAL_NAME, SUPER_PATHWAY) %>% 
  right_join(corr_metab, by = "Var1") -> corr_metab

df_chem %>% 
  mutate(Var2 = CHEM_ID) %>% 
  dplyr::select(Var2, CHEMICAL_NAME, SUPER_PATHWAY) %>% 
  right_join(corr_metab, by = "Var2") -> corr_metab

corr_metab %>% 
  drop_na(p) -> corr_metab


## results available as "correlation_153_testtrain.csv"



###################################
## Step 5: Generate Circos plot ##
##################################

## read association results for different phenotypes of interest
obesity_res <- fread("HELIOS_Metabolon1073_glm_Obesity8129_AgeSexEthnicitySSID_2023.csv", 
                     header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)

ht_res <- fread("HELIOS_Metabolon1073_glm_HT8143_AgeSexEthnicitySSID_2023.csv",
                header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)

t2d_res <- fread("HELIOS_Metabolon1073_glm_T2D8143_AgeSexEthnicitySSID_2023.csv", 
                   header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)

cvd_res <- fread("HELIOS_Metabolon1073_glm_CVD8143_AgeSexEthnicitySSID_2023.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)


## subset the association results for 153 metabolites and add column denoting if significant for the specific phenotype
obesity_res %>% 
  filter(CHEM_ID %in% final_out_sig$CHEM_ID) %>% 
  mutate(sig_obesity = ifelse(adjP_Bonf < 0.05, 0.5, NA)) %>% 
  dplyr::select(CHEM_ID, sig_obesity) -> obesity_res_new

ht_res %>% 
  filter(CHEM_ID %in% final_out_sig$CHEM_ID) %>% 
  mutate(sig_ht = ifelse(adjP_Bonf < 0.05, 0.5, NA)) %>% 
  dplyr::select(CHEM_ID, sig_ht) -> ht_res_new

t2d_res %>% 
  filter(CHEM_ID %in% final_out_sig$CHEM_ID) %>% 
  mutate(sig_t2d = ifelse(adjP_Bonf < 0.05, 0.5, NA)) %>% 
  dplyr::select(CHEM_ID, sig_t2d) -> t2d_res_new

cvd_res %>% 
  filter(CHEM_ID %in% final_out_sig$CHEM_ID) %>% 
  mutate(sig_cvd = ifelse(adjP_Bonf < 0.05, 0.5, NA)) %>% 
  dplyr::select(CHEM_ID, sig_cvd) -> cvd_res_new


## merge results and format into split matrix for circos plot
# merge different results to be plotted
final_out_sig %>% 
  filter(!SUPER_PATHWAY %in% c("Partially Characterized Molecules", "")) %>%
  left_join(ht_res_new, by = "CHEM_ID") %>% 
  left_join(obesity_res_new, by = "CHEM_ID") %>% 
  left_join(t2d_res_new, by = "CHEM_ID") %>% 
  left_join(cvd_res_new, by = "CHEM_ID") %>% 
  arrange(SUPER_PATHWAY, SUB_PATHWAY) %>% 
  dplyr::select(SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, Est_I, Est_M, Est_IM, sig_ht, sig_obesity, sig_t2d, sig_cvd) %>% 
  group_by(SUPER_PATHWAY) %>% 
  mutate(x_index = row_number()) %>% 
  ungroup() -> df_matrix

# save matrix for future use
df_matrix_temp <- df_matrix

# split matrix by super_pathway
split_matrix <- as.factor(df_matrix$SUPER_PATHWAY)

# format as matrix 
df_matrix %>% 
  dplyr::select(-SUPER_PATHWAY, -SUB_PATHWAY) %>% 
  column_to_rownames(var = "CHEMICAL_NAME") -> df_matrix
df_matrix <- as.matrix(df_matrix)


## set plot parameters
circos.clear ()
circos.par(cell.padding = c(0.02, 0, 0.02, 0), start.degree = 90, gap.degree = 1)
circos.par(gap.after = c(1, 1, 1, 1, 1, 1, 7))
set_track_gap(0.005)
track.margin = c(0,0,0,0)


## initialize plot
circos.heatmap.initialize(df_matrix, split = split_matrix, cluster = FALSE)


## create track with super_pathway 
col_SP = c("#e41a1c", "#377eb8", "#ff7f00", "#984ea3", "#4daf4a", "#ffff33", "#f781bf")
circos.track(split_matrix, ylim = c(0, 1), bg.col = col_SP, bg.border = NA,
             panel.fun=function(x, y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$cell.ylim[2] + mm_y(4),
                           CELL_META$sector.index)
             }, track.height = 0.05, cell.padding = c(0.02, 0, 0.02, 0))


## create tracks for disease associations
# HT
circos.track(ylim = c(0, 0.8), bg.border = "grey", bg.col = "#deebf7",
             panel.fun=function(x, y) {
               y = df_matrix[CELL_META$subset, 4]
               x = seq_along(y) - 0.5
               circos.points(x, y,
                             col = "black", pch = 16, cex = 0.5)
             }, track.height = 0.06)


# Obesity
circos.track(ylim = c(0, 0.8), bg.border = "grey", bg.col = "lightyellow", 
             panel.fun=function(x, y) {
               y = df_matrix[CELL_META$subset, 5]
               x = seq_along(y) - 0.5
               circos.points(x, y,
                             col = "black", pch = 16, cex = 0.5)
             }, track.height = 0.06)



# T2D
circos.track(ylim = c(0, 0.8), bg.border = "grey", bg.col = "#deebf7", 
             panel.fun=function(x, y) {
               y = df_matrix[CELL_META$subset, 6]
               x = seq_along(y) - 0.5
               circos.points(x, y,
                             col = "black", pch = 16, cex = 0.5)
             }, track.height = 0.06)


# CVD
circos.track(ylim = c(0, 0.8), bg.border = "grey", bg.col = "lightyellow", 
             panel.fun=function(x, y) {
               y = df_matrix[CELL_META$subset, 7]
               x = seq_along(y) - 0.5
               circos.points(x, y,
                             col = "black", pch = 16, cex = 0.5)
             }, track.height = 0.06)


## create heatmap track for ethnicity association effect sizes
col_est = colorRamp2(c(-1.5, 0, 1.8), c("firebrick3", "white", "blue3"))
circos.heatmap(df_matrix[, 1:3], split = split_matrix, col = col_est, 
               rownames.side = "none", 
               track.height = 0.2,
               bg.border = NA, cell.border = "grey", cell.lwd = 0.2, cell.lty = 1, 
               show.sector.labels = FALSE)


## draw links for correlation between metabolites
#create dataset 
df_matrix_temp %>% 
  mutate(SUPER_PATHWAY.x = SUPER_PATHWAY,
         SUPER_PATHWAY.y = SUPER_PATHWAY,
         SUB_PATHWAY.x = SUB_PATHWAY,
         SUB_PATHWAY.y = SUB_PATHWAY,
         CHEMICAL_NAME.x = CHEMICAL_NAME,
         CHEMICAL_NAME.y = CHEMICAL_NAME, 
         index.x = x_index, 
         index.y = x_index) %>% 
  dplyr::select(ends_with(c(".x", ".y"))) -> temp

temp %>% 
  dplyr::select(ends_with(".x")) %>% 
  right_join(corr_metab, by = c("SUPER_PATHWAY.x", "CHEMICAL_NAME.x")) -> temp_plot

temp %>% 
  dplyr::select(ends_with(".y")) %>% 
  right_join(temp_plot, by = c("SUPER_PATHWAY.y", "CHEMICAL_NAME.y")) -> temp_plot

# plot all correlations with absolute correlation greater than equal to 0.6
temp_plot %>% 
  filter(!SUPER_PATHWAY.x %in% c("Partially Characterized Molecules", "")) %>% 
  filter(!SUPER_PATHWAY.y %in% c("Partially Characterized Molecules", "")) %>% 
  filter(abs(corr) >= 0.6) -> temp_plot
  
for(i in 1:nrow(temp_plot)) {
  circos.link(sector.index1 = temp_plot$SUPER_PATHWAY.x[i], point1 = temp_plot$index.x[i], 
              sector.index2 = temp_plot$SUPER_PATHWAY.y[i], point2 = temp_plot$index.y[i], 
              col = "grey", lwd = 0.8)
}

# plot subset of correlations between metabolites in different super_pathways
temp_plot %>% 
  filter(SUPER_PATHWAY.x != SUPER_PATHWAY.y) -> temp_plot_2

for(i in 1:nrow(temp_plot_2)) {
  circos.link(sector.index1 = temp_plot_2$SUPER_PATHWAY.x[i], point1 = temp_plot_2$index.x[i], 
              sector.index2 = temp_plot_2$SUPER_PATHWAY.y[i], point2 = temp_plot_2$index.y[i], 
              col = "forestgreen", lwd = 0.8)
}

# plot subset of correlations between metabolites in same super_pathway but different sub_pathways
temp_plot %>% 
  filter(SUPER_PATHWAY.x == SUPER_PATHWAY.y & SUB_PATHWAY.x != SUB_PATHWAY.y) -> temp_plot_3

for(i in 1:nrow(temp_plot_3)) {
  circos.link(sector.index1 = temp_plot_3$SUPER_PATHWAY.x[i], point1 = temp_plot_3$index.x[i], 
              sector.index2 = temp_plot_3$SUPER_PATHWAY.y[i], point2 = temp_plot_3$index.y[i], 
              col = "darkblue", lwd = 0.8)
}


## draw legends
lgd_corr = Legend(title = "regression coefficient", col_fun = col_est)
lgd_metab = Legend(title = "metabolite categories", legend_gp = gpar(fill = col_SP), at = unique(as.character(split_matrix)))
lgd_link = Legend(title = "correlated metabolites", legend_gp = gpar(fill = c("grey", "darkblue", "forestgreen")), at = c("same pathway in same category", "different pathway in same category", "different category"))

lgd_list = packLegend(lgd_metab, lgd_corr, lgd_link)
draw(lgd_list, x = unit(0.8, "npc"), y = unit(0.5, "npc"), just = "left")

circos.clear()



###########################################
## Step 6: Metabolic variance explained ##
##########################################

## read dataset with age, sex, bmi, genetic ancestry pcs, and dietary pcs
df_pheno_new <- fread("<insert filename of dataset>", 
                sep = "\t",  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


## subset 153 metabolite data for same number of samples
df_metab_153 %>% 
  filter(PARENT_SAMPLE_NAME %in% df_pheno_new$PARENT_SAMPLE_NAME) %>%
  select_if(~ !any(is.na(.))) -> df_metab_new

# check order of samples
table(df_metab_new$PARENT_SAMPLE_NAME == df_pheno_new$PARENT_SAMPLE_NAME)


## remove sample name columns
df_metab_new %>% 
  column_to_rownames("PARENT_SAMPLE_NAME") -> df_metab_new

df_pheno_new %>% 
  column_to_rownames("PARENT_SAMPLE_NAME") -> df_pheno_new


## partial correlation
merged_res = NULL
blank_df = as.data.frame(t(rep(NA, 6)))
colnames(blank_df) <- c("estimate", "p.value", "statistic", "n", "gp", "Method")

for(i in 1:ncol(df_metab_new)) {
  for(j in 1:ncol(df_pheno_new)) {
    corr_res <- tryCatch({pcor.test(x = df_metab_new[,i], y = df_pheno_new[, j], z = df_pheno_new[, -j], method = "pearson")}, error = function(error) {return(blank_df)})
    variables <- paste(colnames(df_metab_new[i]), colnames(df_pheno_new[j]), sep = " ")
    curr_res <- cbind(corr_res, variables)
    merged_res <- rbind(merged_res, curr_res)
  }
}  


## format results: 
    #get partial rsquared values from correlation estimates 
    #add up rsquared values of 50 genetic pcs for each metabolite 
    #add up rsquared values of 20 dietary pcs  for each metabolites

merged_res %>% 
  separate(variables, into = c("Var1", "Var2"), sep = " ") %>% 
  mutate(rsquared = estimate^2) %>% 
  dplyr::select(Var1, Var2, rsquared) %>% 
  filter(Var2 %in% c("Age", "Sex", "Bmi")) -> res_1

merged_res %>% 
  separate(variables, into = c("Var1", "Var2"), sep = " ") %>% 
  mutate(rsquared = estimate^2) %>% 
  filter(str_detect(Var2, "gPC")) %>% 
  group_by(Var1) %>% 
  summarise(rsquared = sum(rsquared)) %>% 
  mutate(Var2 = "Genetic Ancestry") %>% 
  dplyr::select(Var1, Var2, rsquared) -> res_2

merged_res %>% 
  separate(variables, into = c("Var1", "Var2"), sep = " ") %>% 
  mutate(rsquared = estimate^2) %>% 
  filter(str_detect(Var2, "fPC")) %>% 
  group_by(Var1) %>% 
  summarise(rsquared = sum(rsquared)) %>% 
  mutate(Var2 = "Diet") %>% 
  dplyr::select(Var1, Var2, rsquared) -> res_3

res_ppcor <- rbind(res_1, res_2, res_3)


## plot partial rsquared values for each variable for all metabolites 
res_ppcor %>% 
  mutate(Var2 = ifelse(Var2 == "Bmi", "BMI", Var2)) %>%
  mutate(Var2 = fct_reorder(Var2, rsquared, .fun='mean')) %>%
  ggplot(aes(x = Var2, y = rsquared)) + 
  geom_boxplot(aes(fill = Var2), alpha = 0.8, width = 0.7) + 
  xlab("") + ylab("Partial R-squared") + 
  theme_classic(20) +
  theme(legend.position="none")


