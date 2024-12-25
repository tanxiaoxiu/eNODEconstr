library(tidyverse)
library(ggplot2)
library(patchwork)


############plot FigureS1
N <- 10
M <- 10
compare_path <- "~/eNODEconstr/simulation/comparison/DKI"

path <- paste0(compare_path, "/n", N, "m", M)
setwd(path)
all_results <- readRDS("all_results_DKI.rds")

###########score r_rmse
all_score_r_rmse_mean <- all_results$all_score_r_rmse_mean
all_score_r_rmse_mean <- na.omit(all_score_r_rmse_mean)

all_score_r_rmse_mean$Method <- factor(all_score_r_rmse_mean$Method, levels = c("DKI", "eNODEconstr"))
all_score_r_rmse_mean$SubjectN <- factor(all_score_r_rmse_mean$SubjectN)
all_score_r_rmse_mean$Timepoint <- as.factor(all_score_r_rmse_mean$Timepoint)

subject_labels <- c(
  "10" = "italic(S)~'='~10",
  "15" = "italic(S)~'='~15",
  "20" = "italic(S)~'='~20"
)

plot_score_r_rmse_mean <- ggplot(all_score_r_rmse_mean, aes(x = Timepoint, y = r_rmse, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid( ~ SubjectN, labeller = labeller(SubjectN = as_labeller(subject_labels, label_parsed))) +
  scale_fill_manual(values = c("#82b0d2","#df7676")) +
  labs( y = "Relative RMSE", fill = "Method") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 35,color = "black"),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 33,color = "black", margin = margin(t = 18)), 
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 30,color = "black", margin = margin(l = 18)), 
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    legend.position = "none",  
    legend.title = element_text(size = 34),  
    legend.text = element_text(size = 34)  
  )
ggsave(filename="FigureS1A.png",plot=plot_score_r_rmse_mean,device="png",dpi=600,units="in",width=30,height=7)


###########score spearman
all_score_spearman_mean <- all_results$all_score_spearman_mean
all_score_spearman_mean <- na.omit(all_score_spearman_mean)
all_score_spearman_mean$Method <- factor(all_score_spearman_mean$Method, levels = c("DKI", "eNODEconstr"))
all_score_spearman_mean$SubjectN <- factor(all_score_spearman_mean$SubjectN)
all_score_spearman_mean$Timepoint <- as.factor(all_score_spearman_mean$Timepoint)

plot_score_spearman_mean <- ggplot(all_score_spearman_mean, aes(x = Timepoint, y = spearman, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid( ~ SubjectN) +
  scale_fill_manual(values = c("#82b0d2","#df7676")) +
  labs(x = expression(italic(T)), y = "Spearman correlation", fill = "Method") +
  theme_minimal() +
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(color = "white"),
    axis.title.x = element_text(size = 33,color = "black", margin = margin(t = 18)), 
    axis.title.y = element_text(size = 33,color = "black"),
    axis.text.x = element_text(size = 30,color = "black"),  
    axis.text.y = element_text(size = 30,color = "black", margin = margin(l = 18)), 
    panel.grid.major = element_line(color = "grey90"),
    
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    legend.position = "bottom",  
    legend.title = element_text(size = 33, margin = margin(r = 20)), 
    legend.text = element_text(size = 32, margin = margin(r = 10)),  
    legend.key.size = unit(1, "cm") 
  )

ggsave(filename="FigureS1B.png",plot=plot_score_spearman_mean,device="png",dpi=600,units="in",width=30,height=8)



