library(tidyverse)
library(ggplot2)
library(patchwork)


############plot
N=20
M=20
base_path1 <- "~/eNODEconstr/simulation/"
base_path2 <- paste0(base_path1,"comparison")

path <- paste0(base_path2, "/n", N, "m", M)
setwd(path)
all_results <- readRDS("all_results.rds")


#########traj_r_rmse
all_traj_r_rmse_subject <- all_results$all_traj_r_rmse_subject
all_traj_r_rmse_subject <- na.omit(all_results$all_traj_r_rmse_subject)

all_traj_r_rmse_subject$Method <- factor(all_traj_r_rmse_subject$Method, levels = c("NODE","eNODE","eNODEconstr"))
all_traj_r_rmse_subject$SubjectN <- factor(all_traj_r_rmse_subject$SubjectN)
all_traj_r_rmse_subject$Timepoint <- as.factor(all_traj_r_rmse_subject$Timepoint)
all_traj_r_rmse_subject$Group <- factor(all_traj_r_rmse_subject$Group, levels = c("Microbe", "Metabolite"))

subject_labels <- c(
  "10" = "italic(S)~'='~10",
  "20" = "italic(S)~'='~20",
  "40" = "italic(S)~'='~40"
)


# Create the plot
plot_traj_r_rmse_subject <- ggplot(all_traj_r_rmse_subject, aes(x = Timepoint, y = r_rmse, fill = Method)) +
  geom_boxplot() +
  facet_grid(Group ~ SubjectN, labeller = labeller(SubjectN = as_labeller(subject_labels, label_parsed))) +
  scale_fill_manual(values = c("#3F918B","#896191","#df7676")) +
  labs(x = expression(italic(T)), y = "Relative RMSE", fill = "Method") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 34,color = "black"),
    axis.title.x = element_text(size = 33,color = "black", margin = margin(t = 19)), 
    axis.title.y = element_text(size = 33,color = "black", margin = margin(r = 18)), 
    axis.text.x = element_text(size = 30,color = "black"),  
    axis.text.y = element_text(size = 30,color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    legend.position = "bottom",  
    legend.title = element_text(size = 33, margin = margin(r = 20)), 
    legend.text = element_text(size = 32, margin = margin(r = 10)),  
    legend.key.size = unit(1, "cm")
  )

ggsave(filename="plot_FigureS2a.png",plot=plot_traj_r_rmse_subject,device="png",dpi=600,units="in",width=13,height=15)


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
  "20" = "italic(S)~'='~20",
  "40" = "italic(S)~'='~40"
)

plot_score_r_rmse_subject <- ggplot(all_score_r_rmse_subject, aes(x = Timepoint, y = r_rmse, fill = Method)) +
  geom_boxplot() +
  facet_grid(Group ~ SubjectN, labeller = labeller(SubjectN = as_labeller(subject_labels, label_parsed))) +
  scale_fill_manual(values = c("#3F918B","#896191","#df7676")) +
  labs(x = expression(italic(T)), y = "Relative RMSE", fill = "Method") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 34,color = "black"),
    axis.title.x = element_text(size = 33,color = "black", margin = margin(t = 19)), 
    axis.title.y = element_text(size = 33,color = "black", margin = margin(r = 18)), 
    axis.text.x = element_text(size = 30,color = "black"),  
    axis.text.y = element_text(size = 30,color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    legend.position = "bottom",  
    legend.title = element_text(size = 33, margin = margin(r = 20)), 
    legend.text = element_text(size = 32, margin = margin(r = 10)),  
    legend.key.size = unit(1, "cm")
  )
ggsave(filename="plot_FigureS2b.png",plot=plot_score_r_rmse_subject,device="png",dpi=600,units="in",width=13,height=15)


#########mu
all_mu_index_mean <- readRDS("all_mu_index_mean.rds")
plot_list <- list()
for (subj in subject_values) {
  for (time in timepoint_values) {
    data <- subset(all_mu_index_mean, SubjectN == subj & Timepoint == time)
    data <- data[,c("index","value","Method")]
    dat_wide <- data %>%
      select(index, value, Method) %>%
      pivot_wider(names_from = index, values_from = value)
    
    .plot <- ggradar(
      dat_wide,
      background.circle.transparency = 0, 
      background.circle.colour = "white", 
      group.colours = c("#896191","#df7676"), 
      grid.min = 0,   
      grid.mid = 0.5,  
      grid.max = 1.0, 
      values.radar = c(0, 0.5, 1.0),
      axis.label.size = 10, 
      grid.label.size = 9, 
      gridline.mid.colour = "grey",
      axis.label.offset = 1.15 
    )+
      ggtitle(bquote(italic(S)~"="~.(subj)))+
      labs(color = "Method") +  
      theme(
        plot.title = element_text(hjust = 0.5, size = 28),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),                
        legend.position = "bottom"                            
      )
    
    plot_list[[paste("Subject", subj, "Timepoint", time)]] <- .plot
  }
}


combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 1, guides = "collect") &
  theme(
    legend.position = "bottom",             
    legend.text = element_text(size = 29),  
    legend.key.size = unit(4, "cm"),        
    legend.title = element_text(size = 30),
    plot.margin = unit(c(0.1, -0.5, 0.1, -0.5), "cm")
  )

ggsave(filename="plot_FigureS2c.png",plot=combined_plot,device="png",dpi=600,units="in",width=24,height=8)

