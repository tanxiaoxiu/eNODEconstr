library(tidyverse)
library(ggplot2)
library(patchwork)


############plot FigureS3
N=10
M=10
base_path1 <- "~/eNODEconstr/simulation/"
base_path2 <- paste0(base_path1,"comparison")
path <- paste0(base_path2, "/n", N, "m", M)
setwd(path)
all_results <- readRDS("all_results_t3.rds")

#s########traj_r_rmse
all_traj_r_rmse_subject <- all_results$all_traj_r_rmse_subject
all_traj_r_rmse_subject <- na.omit(all_results$all_traj_r_rmse_subject)

all_traj_r_rmse_subject$Method <- factor(all_traj_r_rmse_subject$Method, levels = c("NODE","eNODE","eNODEconstr"))
all_traj_r_rmse_subject$SubjectN <- factor(all_traj_r_rmse_subject$SubjectN)
all_traj_r_rmse_subject$Timepoint <- as.factor(all_traj_r_rmse_subject$Timepoint)
all_traj_r_rmse_subject$Group <- factor(all_traj_r_rmse_subject$Group, levels = c("Microbe", "Metabolite"))

subject_labels <- c(
  "10" = "italic(S)~'='~10",
  "15" = "italic(S)~'='~15",
  "20" = "italic(S)~'='~20",
  "40" = "italic(S)~'='~40",
  "60" = "italic(S)~'='~60"
)

plot_traj_r_rmse_subject <- ggplot(all_traj_r_rmse_subject, aes(x = Timepoint, y = r_rmse, fill = Method)) +
  geom_boxplot() +
  facet_grid(Group ~ SubjectN, labeller = labeller(SubjectN = as_labeller(subject_labels, label_parsed))) +
  scale_fill_manual(values = c("#13679e","#f09b81","#df7676")) +
  labs(x = expression(italic(T)), y = "Relative RMSE", fill = "Method") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 35,color = "black"),
    axis.title.x = element_text(size = 33,color = "black"), 
    axis.title.y = element_text(size = 33,color = "black"),
    axis.text.x = element_text(size = 30,color = "black"),  
    axis.text.y = element_text(size = 30,color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    legend.position = "bottom",  
    legend.title = element_text(size = 34),  
    legend.text = element_text(size = 34)  
  )
ggsave(filename="FigureS3a.png",plot=plot_traj_r_rmse_subject,device="png",dpi=600,units="in",width=18,height=15)


#########Score_r_rmse
all_score_r_rmse_subject <- all_results$all_score_r_rmse_subject
all_score_r_rmse_subject <- na.omit(all_results$all_score_r_rmse_subject)

all_score_r_rmse_subject$Method <- factor(all_score_r_rmse_subject$Method, levels = c("NODE","eNODE","eNODEconstr"))
all_score_r_rmse_subject$SubjectN <- factor(all_score_r_rmse_subject$SubjectN)
all_score_r_rmse_subject$Timepoint <- as.factor(all_score_r_rmse_subject$Timepoint)
all_score_r_rmse_subject$Group <- gsub("^Microbe$", "Composition score", all_score_r_rmse_subject$Group)
all_score_r_rmse_subject$Group <- gsub("^Metabolite$", "Function score", all_score_r_rmse_subject$Group)
all_score_r_rmse_subject$Group <- factor(all_score_r_rmse_subject$Group, levels = c("Composition score", "Function score"))

subject_labels <- c(
  "10" = "italic(S)~'='~10",
  "15" = "italic(S)~'='~15",
  "20" = "italic(S)~'='~20",
  "40" = "italic(S)~'='~40",
  "60" = "italic(S)~'='~60"
)

plot_score_r_rmse_subject <- ggplot(all_score_r_rmse_subject, aes(x = Timepoint, y = r_rmse, fill = Method)) +
  geom_boxplot() +
  facet_grid(Group ~ SubjectN, labeller = labeller(SubjectN = as_labeller(subject_labels, label_parsed))) +
  scale_fill_manual(values = c("#13679e","#f09b81","#df7676")) +
  labs(x = expression(italic(T)), y = "Relative RMSE", fill = "Method") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 35,color = "black"),
    axis.title.x = element_text(size = 33,color = "black"), 
    axis.title.y = element_text(size = 33,color = "black"),
    axis.text.x = element_text(size = 30,color = "black"),  
    axis.text.y = element_text(size = 30,color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    legend.position = "bottom",  
    legend.title = element_text(size = 34),
    legend.text = element_text(size = 34) 
  )
ggsave(filename="FigureS3b.png",plot=plot_score_r_rmse_subject,device="png",dpi=600,units="in",width=18,height=15)


#########mu
all_mu_index_mean <- readRDS("all_mu_index_mean_t3.rds")
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
      group.colours = c("#f09b81","#df7676"), 
      grid.min = 0,    
      grid.mid = 0.5,  
      grid.max = 1.0, 
      values.radar = c(0, 0.5, 1.0),
      axis.label.size = 8, 
      grid.label.size = 7, 
      gridline.mid.colour = "grey",
      axis.label.offset = 1.15 
    )+
      ggtitle(bquote(italic(S)~"="~.(subj)))+
      
      theme(
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 28)
        
      )
    
    plot_list[[paste("Subject", subj, "Timepoint", time)]] <- plot
  }
}

combined_plot <- wrap_plots(plot_list, ncol = 5, nrow = 1) &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 26)
  ) 

combined_plot_with_legend <- combined_plot +
  plot_layout(guides = "collect")&
  theme(
    plot.margin = unit(c(0.2, -0.5, 0.2, -0.5), "cm")
  )
combined_plot_with_legend
ggsave(filename="FigureS3c.png",plot=combined_plot_with_legend,device="png",dpi=600,units="in",width=28,height=6)



