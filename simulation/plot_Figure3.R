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


############plot Figure3
N=10
M=10
s_values <- c(10,15,20)
t_values <- c(5, 10, 15, 20, 25)
base_folders <- c("NODE","eNODE","eNODEconstr")
Iter <- 10
path <- paste0(base_path2, "/n", N, "m", M)
setwd(path)
all_results <- readRDS("all_results.rds")


base_path1 <- "~/eNODEconstr/simulation/"
base_path2 <- paste0(base_path1,"comparison")

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
  "20" = "italic(S)~'='~20"
)

# Create the plot
plot_score_r_rmse_subject <- ggplot(all_score_r_rmse_subject, aes(x = Timepoint, y = r_rmse, fill = Method)) +
  geom_boxplot() +
  facet_grid(Group ~ SubjectN, labeller = labeller(SubjectN = as_labeller(subject_labels, label_parsed))) +
  scale_fill_manual(values = c("#3F918B","#896191","#df7676")) +
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

ggsave(filename="Figure3.png",plot=plot_score_r_rmse_subject,device="png",dpi=600,units="in",width=30,height=15)


