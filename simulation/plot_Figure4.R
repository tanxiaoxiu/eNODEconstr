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
N=10
M=10
Savepath <- paste0(base_path2, "/n", N, "m", M)

setwd(Savepath)

base_folders <- c("eNODE","eNODEconstr")
subject_values <- c(10, 15, 20)
timepoint_values <- c(5, 10, 15, 20, 25)
all_mu_index_mean <- readRDS("all_mu_index_mean.rds")

plot_list <- list()
plot_list <- list()
for (subj in subject_values) {
  for (time in timepoint_values) {
    data <- subset(all_mu_index_mean, SubjectN == subj & Timepoint == time)
    data <- data[,c("index", "value", "Method")]
    dat_wide <- data %>%
      select(index, value, Method) %>%
      pivot_wider(names_from = index, values_from = value)
    
    .plot <- ggradar(
      dat_wide,
      background.circle.transparency = 0, 
      background.circle.colour = "white", 
      group.colours = c("#896191", "#df7676"), 
      grid.min = 0,    
      grid.mid = 0.5,  
      grid.max = 1.0, 
      values.radar = c(0, 0.5, 1.0),
      axis.label.size = 7.5, 
      grid.label.size = 6.5, 
      gridline.mid.colour = "grey",
      axis.label.offset = 1.15 
    ) +
      labs(color = "Method") +  
      theme(
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),                
        legend.position = "bottom"                           
      )
    
    plot_list[[paste("Subject", subj, "Timepoint", time)]] <- .plot
  }
}

combined_plot <- wrap_plots(plot_list, ncol = 5, nrow = 3, guides = "collect") &
  theme(
    legend.position = "bottom",             
    legend.text = element_text(size = 28),  
    legend.key.size = unit(4, "cm"),        
    legend.title = element_text(size = 30) 
  )


ggsave(filename="Figure4.png",plot=combined_plot,device="png",dpi=600,units="in",width=25,height=13)



