

### This script consolidate the code use to plot the enrichment plot and heatmap

#1. PCA heatmap

library(ggplot2)
library(reshape2)
library(stringr)

#Read the file that contain the beta effect and p value to generate heatmap
PCA <- read.csv("Epigenetic_PC1_5_regression_beta_summary_trimed_methods_heatmap_p05.csv")

# Step 1: Subset the data frame to include only the columns of interest
subset_df <- PCA[, c("Group","Traits","PC1_Estimate", "PC2_Estimate", "PC3_Estimate", "PC4_Estimate", "PC5_Estimate",
                    "PC1_P", "PC2_P", "PC3_P", "PC4_P", "PC5_P")]

melt_df <- melt(subset_df, id = c("Group", "Traits"))

melt_df_estimate <- melt_df %>% filter(variable %in% c("PC1_Estimate", "PC2_Estimate", "PC3_Estimate", "PC4_Estimate", "PC5_Estimate"))
melt_df_p <- melt_df %>% filter(variable %in% c("PC1_P", "PC2_P", "PC3_P", "PC4_P", "PC5_P"))

melt_df_final <- as.data.frame(cbind(melt_df_estimate,melt_df_p))

names(melt_df_final) <- c("Group","Variable","variable","b_estimate","Group1","Variable1","variable1","p_value")

#Assign * and ** for p <0.05 and p < 0.0015
melt_df_final$signif <- ifelse(melt_df_final$p_value>0.05,NA,
                               ifelse(melt_df_final$p_value<=0.0015,"**","*"))

#filter the data into two drop to plot two heatmap
melt_df_final_plot1 <- melt_df_final %>% filter(Group %in% c("Anthropometry","Lung capacity","Lab measures"))
melt_df_final_plot2 <- melt_df_final %>% filter(Group %in% c("Cognition","Diet: Dietary score & metabolite score","Diet: Perkcal food group","Cardiovascular","PRS"))

#plot first heatmap
png("~/Desktop/tmp1.png",height = 12,width = 50,units = "cm",res = 360)
ggplot(data = melt_df_final_plot1, aes(x = Variable, y = variable, fill = b_estimate)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "#FF0000",
                       limits = c(-0.5, 0.5),
                       guide = guide_colorbar(title.position = "top", # Adjust position and title as needed
                                              title.vjust = 1,        # Adjust vertical position of title
                                              title.theme = element_text(size = 16), # Adjust title font size
                                              label.theme = element_text(size = 14))) + # Adjust label font size) +
  facet_grid(.~Group,scales = "free",space = "free", labeller = label_wrap_gen(width = 2, multi_line = TRUE))+
  #facet_grid(.~Group,scales = "free",space = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        strip.text = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  labs(title = "",x = "", y = "Estimate (PCs)")+
  #scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = signif),size = 6, color = "black", vjust = 0.5, hjust = 0.5) + 
  scale_y_discrete(limits = c("PC5_Estimate","PC4_Estimate","PC3_Estimate","PC2_Estimate","PC1_Estimate"))  # Flip the x-axis

dev.off()

#plot 2nd heatmap
png("~/Desktop/tmp2.png",height = 14.5,width = 50,units = "cm",res = 360)
ggplot(data = melt_df_final_plot2, aes(x = Variable, y = variable, fill = b_estimate)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "#FF0000",
                       limits = c(-0.5, 0.5),
                       guide = guide_colorbar(title.position = "top", # Adjust position and title as needed
                                              title.vjust = 1,        # Adjust vertical position of title
                                              title.theme = element_text(size = 16), # Adjust title font size
                                              label.theme = element_text(size = 14))) + # Adjust label font size) +
  #facet_grid(.~Group,scales = "free",space = "free", labeller = label_wrap_gen(width = 2, multi_line = TRUE))+
  facet_grid(.~Group,scales = "free",space = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        strip.text = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  labs(title = "",x = "", y = "Estimate (PCs)")+
  #scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = signif),size = 6, color = "black", vjust = 0.5, hjust = 0.5) + 
  scale_y_discrete(limits = c("PC5_Estimate","PC4_Estimate","PC3_Estimate","PC2_Estimate","PC1_Estimate"))  # Flip the x-axis

dev.off()

#2. Enrichment plot

library(dplyr)
library(ggplot2)
library(stringr)
library(ggstar)

#read the file that contains the enrichment fold as well as p value
a3 <- read.csv("EWAS_enrichment_Method_paper.csv")

#plot the enrichment plot for 8 selected categories
png("EWAS_enrichment.png", height = 20, width = 30, units = "cm", res = 480)  # Open PNG device for plotting
ggplot(a3 %>% 
         group_by(grp) %>% 
         arrange(fold_enrich) %>%
         mutate(variable = factor(variable, levels = variable)),
       aes(x = fold_enrich, y = variable)) +
  geom_bar(stat = "identity", width = 0.2, aes(fill = grp), color = "white") +  # Add bars
  geom_point(aes(x = fold_enrich, y = variable, color = grp), shape=21, fill="white", position="stack",stat="identity", size=2, show.legend = F)+
  #geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +  # Add vertical line at x = 1
  theme_bw() +  # Apply black and white theme
  facet_wrap(~ grp, ncol = 2, scales = "free_y") + 
  labs(x="",y="Phenotypes")+
  theme(
    legend.position = "none",  # Remove legend
    strip.text.x = element_blank(),  # Remove x-axis strip text
    strip.background = element_rect(fill = "white", colour = "white"),  # White strip background
    panel.spacing = unit(1, "cm"),  # Increase spacing between facets
    axis.text = element_text(size = 16),  # Set text size for axes
    axis.title.x = element_text(size = 20, face = "bold"),  # Set x-axis title size and bold
    axis.title.y = element_blank(),  # Remove y-axis title
    legend.background = element_rect(linewidth = 0.5, colour = "black"),  # Set legend background
    legend.title = element_text(size = 20, face = "bold"),  # Set legend title size and bold
    legend.text = element_text(size = 16),  # Set legend text size
    strip.text.y = element_text(margin = margin(0, 12, 0, ifelse(c(1, 3, 5) %in% unique(a3$grp), 6, 3))) # Adjust margin for specific facets
  ) + 
  coord_cartesian(xlim = c(1, 3.2))

dev.off()  # Close the PNG device

#read the file that contains the enrichment fold as well as p value for PRS
a4 <- read.csv("Enrichment_16K_summary_PRS_sig_Method_paper.csv")

#a4$Signif <- ifelse(a4$Phyper_16K>0.05,"non-sig",
#                    ifelse(a4$Phyper_16K<0.001,"p < 0.001", "p < 0.05"))
#a4$Signif <- factor(a4$Signif,levels = c("non-sig", "p < 0.05", "p < 0.001"))

#start of new plot

list <- as.vector(unique(a4$grp))
#list <- list[c(9,6,4,8,11,5,7,12,10,1,13,3,2)]
a4 %>% group_by(grp) %>% summarise(count=n())
a4$grp <- factor(a4$grp,levels = list)

#subset the results into two group to plot the two enrichment plots
a4_grp1 <- a4 %>% subset(group == 1)
a4_grp2 <- a4 %>% subset(group == 2)

png("EWAS_enrichment_PRS_grp1.png", height = 10, width = 15, units = "cm", res = 480)
ggplot(a4_grp1 %>% 
         group_by(grp) %>% 
         arrange(fold_enrich) %>%
         mutate(variable = factor(variable, levels = variable)),
       aes(x = fold_enrich, y = variable)) +
  geom_bar(stat = "identity", width = 0.3, aes(fill = factor(grp, levels = c("PRS: Body composition", "PRS: Cardiovascular status", "PRS: Cognition", "PRS: Lab measurement", "PRS: Lung capacity", "PRS: Public health"))), color = "white") +  # Add bars
  geom_point(aes(color = factor(grp, levels = c("PRS: Body composition", "PRS: Cardiovascular status", "PRS: Cognition", "PRS: Lab measurement", "PRS: Lung capacity", "PRS: Public health"))), shape = 21, fill = "white", position = position_dodge(width = 0.2), size = 2) +  # Add points
  theme_bw() +
  scale_fill_manual(name = "Type of PRS", values = c("PRS: Body composition" = "#d73027", "PRS: Cardiovascular status" = "#91bfdb", "PRS: Cognition" = "#fee090", "PRS: Lab measurement" = "#7fbf7b", "PRS: Lung capacity" = "#fc8d59", "PRS: Public health" = "#4575b4", guide = guide_legend(reverse=TRUE))) +  # Specify colors for groups
  scale_color_manual(name = "Type of PRS", values = c("PRS: Body composition" = "#d73027", "PRS: Cardiovascular status" = "#91bfdb", "PRS: Cognition" = "#fee090", "PRS: Lab measurement" = "#7fbf7b", "PRS: Lung capacity" = "#fc8d59", "PRS: Public health" = "#4575b4", guide = guide_legend(reverse=TRUE))) +  # Specify colors for groups
  theme(
    strip.text.x = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white"),
    axis.text = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",  # Remove legend
    legend.background = element_blank(),  # Remove legend background
    plot.margin = margin(0.1, 0.5, 0.1, 0.1, "cm")  # Add 1-inch border
    ) + 
  coord_cartesian(xlim = c(1, 2))
dev.off()

png("EWAS_enrichment_PRS_grp2.png", height = 10, width = 15, units = "cm", res = 480)
ggplot(a4_grp2 %>% 
         group_by(grp) %>% 
         arrange(fold_enrich) %>%
         mutate(variable = factor(variable, levels = variable)),
       aes(x = fold_enrich, y = variable)) +
  geom_bar(stat = "identity", width = 0.3, aes(fill = factor(grp, levels = c("PRS: Body composition", "PRS: Cardiovascular status", "PRS: Cognition", "PRS: Lab measurement", "PRS: Lung capacity", "PRS: Public health"))), color = "white") +  # Add bars
  geom_point(aes(color = factor(grp, levels = c("PRS: Body composition", "PRS: Cardiovascular status", "PRS: Cognition", "PRS: Lab measurement", "PRS: Lung capacity", "PRS: Public health"))), shape = 21, fill = "white", position = position_dodge(width = 0.2), size = 2) +  # Add points
  theme_bw() +
  scale_fill_manual(name = "Type of PRS", values = c("PRS: Body composition" = "#d73027", "PRS: Cardiovascular status" = "#91bfdb", "PRS: Cognition" = "#fee090", "PRS: Lab measurement" = "#7fbf7b", "PRS: Lung capacity" = "#fc8d59", "PRS: Public health" = "#4575b4", guide = guide_legend(reverse=TRUE))) +  # Specify colors for groups
  scale_color_manual(name = "Type of PRS", values = c("PRS: Body composition" = "#d73027", "PRS: Cardiovascular status" = "#91bfdb", "PRS: Cognition" = "#fee090", "PRS: Lab measurement" = "#7fbf7b", "PRS: Lung capacity" = "#fc8d59", "PRS: Public health" = "#4575b4", guide = guide_legend(reverse=TRUE))) +  # Specify colors for groups
  theme(
    strip.text.x = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white"),
    axis.text = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",  # Remove legend
    legend.background = element_blank(),  # Remove legend background
    plot.margin = margin(0.1, 0.5, 0.1, 0.1, "cm")  # Add 1-inch border
    ) + 
  coord_cartesian(xlim = c(1, 2))
dev.off()

