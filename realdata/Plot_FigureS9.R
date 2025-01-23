library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)


############FigureS9
stewd("~/eNODEconstr/realdata/Fiber/results/trajectory") 
filter_fiber_SCFA <- read.delim("filter_fiber_SCFA.txt", row.names = 1, sep = '\t', check.names = FALSE)
df_microbe <- filter_fiber_diet[, 5:(N1+N2+5)]

#FigureS9a
Con_microbe <- df_microbe %>% filter(Group == "Control")
Con_microbe$Sample <- rownames(Con_microbe)
Con_microbe <- Con_microbe[,-1]
df_long <- pivot_longer(Con_microbe, cols = -Sample, names_to = "Microbe", values_to = "Abundance")

df_long <- df_long %>%
  group_by(Microbe) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  left_join(df_long, by = "Microbe") %>%
  mutate(Microbe = factor(Microbe, levels = unique(Microbe))) 

plot_Con_microbe <- ggplot(df_long, aes(x = Microbe, y = Abundance)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  ylab("Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = "grey80", size = 0.5))

ggsave(filename="FigureS9a.png",plot=plot_Con_microbe,device="png",dpi=600,units="in",width=20,height=10)


#FigureS9b
Rs_microbe <- df_microbe %>% filter(Group == "Resistant starch")
Rs_microbe$Sample <- rownames(Rs_microbe)
Rs_microbe <- Rs_microbe[,-1]

df_long <- pivot_longer(Rs_microbe, cols = -Sample, names_to = "Microbe", values_to = "Abundance")
df_long <- df_long %>%
  group_by(Microbe) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  left_join(df_long, by = "Microbe") %>%
  mutate(Microbe = factor(Microbe, levels = unique(Microbe))) 

plot_Rs_microbe <- ggplot(df_long, aes(x = Microbe, y = Abundance)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  ylab("Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = "grey80", size = 0.5))

ggsave(filename="FigureS9b.png",plot=plot_Rs_microbe,device="png",dpi=600,units="in",width=20,height=10)


#FigureS9c
In_microbe <- df_microbe %>% filter(Group == "Inulin")
In_microbe$Sample <- rownames(In_microbe)
In_microbe <- In_microbe[,-1]
df_long <- pivot_longer(In_microbe, cols = -Sample, names_to = "Microbe", values_to = "Abundance")

df_long <- df_long %>%
  group_by(Microbe) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  left_join(df_long, by = "Microbe") %>%
  mutate(Microbe = factor(Microbe, levels = unique(Microbe))) 

plot_In_microbe <- ggplot(df_long, aes(x = Microbe, y = Abundance)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  ylab("Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = "grey80", size = 0.5))

ggsave(filename="FigureS9c.png",plot=plot_In_microbe,device="png",dpi=600,units="in",width=20,height=10)


#FigureS9d
N <- 20
M <- 6
base_path <- "~/eNODEconstr/realdata/Fiber/results/trajectory"
true_path <- "~/eNODEconstr/realdata/Fiber/data"
setwd(true_path)
abundance_data <- read.table("filter_fiber_SCFA.txt",header = T,sep="\t")
microbe_name <- as.data.frame(colnames(abundance_data)[6:(N+5)])
colnames(microbe_name) <- "Microbe_name"
microbe_name$Regulate <- paste0("Microbe",rep(1:N))

#Control
path <- paste0(base_path, "/Control")
setwd(path)
Score <- read.table("score_normalized_pre_subject_l2.txt",header = T,sep="\t")
Score <- Score %>%
  left_join(microbe_name, by = c("Regulate" = "Regulate")) %>%
  mutate(Regulate = Microbe_name) %>%
  select(-Microbe_name)

#Microbe
microbe_data  <- subset(Score, Group == "Microbe")
metabolite_data  <- subset(Score, Group == "Metabolite")
result_df <- data.frame(Subject = character(), SpearmanCorrelation = numeric(), stringsAsFactors = FALSE)
subjects <- unique(microbe_data$Subject)

for (subject in subjects) {
  microbe_sub <- subset(microbe_data, Subject == subject)
  metabolite_sub <- subset(metabolite_data, Subject == subject)
  microbe_sub <- microbe_sub[order(microbe_sub$Regulate, microbe_sub$Score), ]
  metabolite_sub <- metabolite_sub[order(metabolite_sub$Regulate, metabolite_sub$Score), ]
  correlation <- cor(microbe_sub$Score, metabolite_sub$Score, method = "spearman")
  result_df <- rbind(result_df, data.frame(Subject = subject, SpearmanCorrelation = correlation))
}

p1 <-ggplot(result_df, aes(x = "Con", y = SpearmanCorrelation)) +
  geom_boxplot(fill = "#87CEEB", color = "darkblue", outlier.shape = NA, width = 0.3) +  
  geom_jitter(aes(color = SpearmanCorrelation), width = 0.2, size = 3, shape = 16, alpha = 0.7) +  
  scale_color_gradient(low = "blue", high = "red") + 
  ylab("Spearman Correlation") +
  theme_minimal(base_size = 15) +  
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.title.x = element_blank(),  
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25,  color = "black"),
        legend.position = "none",  
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = "grey80", size = 0.5))
#ggsave(filename="Spearman_Microbe_Metabolite_Scores_Subject_Con.png",plot=p1,device="png",dpi=600,units="in",width=8,height=5)


#Rs
path <- paste0(base_path, "/Rs")
setwd(path)
Score <- read.table("score_normalized_pre_subject_l2.txt",header = T,sep="\t")
Score <- Score %>%
  left_join(microbe_name, by = c("Regulate" = "Regulate")) %>%
  mutate(Regulate = Microbe_name) %>%
  select(-Microbe_name)

#Microbe
microbe_data  <- subset(Score, Group == "Microbe")
metabolite_data  <- subset(Score, Group == "Metabolite")
result_df <- data.frame(Subject = character(), SpearmanCorrelation = numeric(), stringsAsFactors = FALSE)
subjects <- unique(microbe_data$Subject)

for (subject in subjects) {
  microbe_sub <- subset(microbe_data, Subject == subject)
  metabolite_sub <- subset(metabolite_data, Subject == subject)
  microbe_sub <- microbe_sub[order(microbe_sub$Regulate, microbe_sub$Score), ]
  metabolite_sub <- metabolite_sub[order(metabolite_sub$Regulate, metabolite_sub$Score), ]
  correlation <- cor(microbe_sub$Score, metabolite_sub$Score, method = "spearman")
  result_df <- rbind(result_df, data.frame(Subject = subject, SpearmanCorrelation = correlation))
}

p2 <-ggplot(result_df, aes(x = "Rs", y = SpearmanCorrelation)) +
  geom_boxplot(fill = "#87CEEB", color = "darkblue", outlier.shape = NA, width = 0.3) +  
  geom_jitter(aes(color = SpearmanCorrelation), width = 0.2, size = 3, shape = 16, alpha = 0.7) +  
  scale_color_gradient(low = "blue", high = "red") + 
  ylab("Spearman Correlation") +
  theme_minimal(base_size = 15) +  
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.title.x = element_blank(),  
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25,  color = "black"),
        legend.position = "none",  
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = "grey80", size = 0.5))
#ggsave(filename="Spearman_Microbe_Metabolite_Scores_Subject_RS.png",plot=p2,device="png",dpi=600,units="in",width=8,height=5)

#In
path <- paste0(base_path, "/In")
setwd(path)
Score <- read.table("score_normalized_pre_subject_l2.txt",header = T,sep="\t")
Score <- Score %>%
  left_join(microbe_name, by = c("Regulate" = "Regulate")) %>%
  mutate(Regulate = Microbe_name) %>%
  select(-Microbe_name)

#Microbe
microbe_data  <- subset(Score, Group == "Microbe")
metabolite_data  <- subset(Score, Group == "Metabolite")
result_df <- data.frame(Subject = character(), SpearmanCorrelation = numeric(), stringsAsFactors = FALSE)
subjects <- unique(microbe_data$Subject)

for (subject in subjects) {
  microbe_sub <- subset(microbe_data, Subject == subject)
  metabolite_sub <- subset(metabolite_data, Subject == subject)
  microbe_sub <- microbe_sub[order(microbe_sub$Regulate, microbe_sub$Score), ]
  metabolite_sub <- metabolite_sub[order(metabolite_sub$Regulate, metabolite_sub$Score), ]
  correlation <- cor(microbe_sub$Score, metabolite_sub$Score, method = "spearman")
  result_df <- rbind(result_df, data.frame(Subject = subject, SpearmanCorrelation = correlation))
}

p3 <-ggplot(result_df, aes(x = "In", y = SpearmanCorrelation)) +
  geom_boxplot(fill = "#87CEEB", color = "darkblue", outlier.shape = NA, width = 0.3) +  
  geom_jitter(aes(color = SpearmanCorrelation), width = 0.2, size = 3, shape = 16, alpha = 0.7) +  
  scale_color_gradient(low = "blue", high = "red") + 
  ylab("Spearman Correlation") +
  theme_minimal(base_size = 15) +  
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.title.x = element_blank(),  
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25,  color = "black"),
        legend.position = "none",  
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = "grey80", size = 0.5))
#ggsave(filename="Spearman_Microbe_Metabolite_Scores_Subject_In.png",plot=p3,device="png",dpi=600,units="in",width=8,height=5)
combined_plot <- grid.arrange(p1, p2, p3, ncol = 3)
setwd(base_path)
ggsave(filename="FigureS9d.png", plot=combined_plot, device="png", dpi=600, units="in", width=12, height=5)




