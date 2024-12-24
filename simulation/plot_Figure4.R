library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(igraph)
library(reticulate)
library(tidyr)
library(cowplot)
library(purrr)
library(data.table)
library(patchwork)
library(ggradar)

#########Plot Figure4
base_path1 <- "~/eNODEconstr/simulation/"
base_path2 <- paste0(base_path1,"comparison")
Savepath <- paste0(base_path2, "/n", N, "m", M)

setwd(Savepath)

base_folders <- c("eNODE","eNODEconstr")
subject_values <- c(10, 15, 20)
timepoint_values <- c(5, 10, 15, 20, 25)
all_mu_index_mean <- readRDS("all_mu_index_mean.rds")

plot_list <- list()
for (subj in subject_values) {
  for (time in timepoint_values) {
    data <- subset(all_mu_index_mean, SubjectN == subj & Timepoint == time)
    data <- data[,c("index","value","Method")]
    dat_wide <- data %>%
      select(index, value, Method) %>%
      pivot_wider(names_from = index, values_from = value)
    
    plot <- ggradar(
      dat_wide,
      background.circle.transparency = 0, 
      background.circle.colour = "white", 
      group.colours = c("#896191","#df7676"), 
      grid.min = 0,    
      grid.mid = 0.5,  
      grid.max = 1.0, 
      values.radar = c(0, 0.5, 1.0),
      axis.label.size = 6, 
      grid.label.size = 5, 
      gridline.mid.colour = "grey",
      axis.label.offset = 1.15 
    )+
      theme(
        legend.position = "none" 
      )
    
    plot_list[[paste("Subject", subj, "Timepoint", time)]] <- plot
  }
}

combined_plot <- wrap_plots(plot_list, ncol = 5, nrow = 3) &
  theme(legend.position = "bottom") 
        

combined_plot_with_legend <- combined_plot +
  plot_layout(guides = "collect")&
  theme(
    plot.margin = unit(c(0.3, 0, 0.3, 0), "cm")) 
  

ggsave(filename="Figure4.png",plot=combined_plot_with_legend,device="png",dpi=600,units="in",width=25,height=15)



