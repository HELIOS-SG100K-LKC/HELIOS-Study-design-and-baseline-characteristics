### This script consolidate the code use to plot the scatterplot, heatmap and enrichment plot

## 1. PCA scatterplot with PRS scatterplot
## 2. Heatmap
## 3. Enrichment plot

## 1. plot PCA scatter plot with beta estimate from PRS as examples

library(ggplot2)
library(dplyr)
library(cowplot)
library(magick)

# Read the PRS beta estimate results
data1 <- read.csv("sentrixid_PRS_PCA.csv")
data1$Ethnicity <- as.factor(data1$Ethnicity)

#filter out participants that are others
data1_new <- data1 %>% filter(Ethnicity != "Others")

#read the PCs regression results
data2 <- read.csv("PC1_PC2_regression_beta_summary_new_trimed_p0.0015.csv") 

#subset them into different group for plotting
data2_clinical <- data2 %>% subset(Group == "Clinical traits")
data2_metabolon <- data2 %>% subset(Group == "Metabolite score") 
data2_PRS <- data2 %>% subset(Group == "PRS")

# Create the first scatter plot (data1) that contains the PCs
p1 <- ggplot(data1_new, aes(x = PC1, y = PC2, color = Ethnicity)) +
  geom_point(aes(fill = Ethnicity), shape = 21, size = 3, stroke = 0.2) +  # Specify fill color based on Ethnicity
  scale_fill_manual(values = c("#fc8d62", "#8da0cb","#66c2a5"),
                    labels = c("Chinese", "Indian", "Malay")) +  # Specify fill colors for legend
  scale_color_manual(values = c("black", "black", "black"),
                     labels = c("Chinese", "Malay", "Indian")) + 
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = "none", # remove legend
    axis.title.x = element_blank(), # remove x-axis label
    axis.title.y = element_blank(), # remove y-axis labelaxis.line = element_line(), # show x and y axes
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"), # customize axis ticks
    # Set the dimensions of the plot
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), # Remove default margins
    plot.title = element_text(size = rel(0)), # Remove default title
    legend.key = element_blank()
    ) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))
ggsave(filename = 'plot1_PCA_blackoutline.png', device = 'png', bg = 'transparent', width = 5, height = 5)

# Create the second scatter plot (data2) contain the PRS results
p2_PRS <- ggplot(data2_PRS, aes(x = Beta_PC1, y = Beta_PC2)) +
  #geom_point() +
  geom_segment(aes(x = 0, y = 0, xend = Beta_PC1, yend = Beta_PC2), color = "red", arrow = arrow(length = unit(0.3, "cm")), linewidth = 1.5) + # Add red arrows
  geom_text_repel(aes(label = Variable_ID), force = 8, size = 4, nudge_x = -0.1, nudge_y = 0.05, fontface = "bold", segment.alpha = 0) + 
  scale_x_continuous(position = "top", limits = c(-0.55, 0.55)) +  # Set x-axis limits and position
  scale_y_continuous(position = "right", limits = c(-0.55, 0.55)) +  # Set y-axis limits and position
  theme(
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plotpanel.grid.major = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = "none",  # Remove legend
    axis.title.x = element_blank(), # remove x-axis label
    axis.title.y = element_blank(), # remove y-axis labelaxis.line = element_line(), # show x and y axes
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"), # customize axis ticks
    # Set the dimensions of the plot
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), # Remove default margins
    plot.title = element_text(size = rel(0)), # Remove default title
  ) +
  coord_cartesian(xlim = c(-0.55, 0.55), ylim = c(-0.55, 0.55))
ggsave(filename = 'plot2_PRS.png', device = 'png', bg = 'transparent', width = 5, height = 5)


#2. PCA heatmap

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

#3. Enrichment plot

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

