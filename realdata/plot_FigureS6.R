library(dplyr)
library(reshape2)
library(ggplot2)


##########FigureS6
#FigureS6a
setwd("~/eNODEconstr/realdata/Synthetic/results/trajectory")
gama_df <- read.table("gama.txt",header = F,sep="\t")
gama_df[1,] <- c("gama","1e0","1e-1", "1e-2","1e-3","1e-4","1e-5","1e-6")
gama_df <- as.data.frame(t(gama_df)) 
colnames(gama_df) <- gama_df[1,]
gama_df <- gama_df[-1,]
gama_df$gama <- factor(gama_df$gama, levels = c("1e0","1e-1", "1e-2","1e-3","1e-4","1e-5","1e-6"))
gama_df$value <- as.numeric(gama_df$value)

p1 <- ggplot(gama_df, aes(x = gama, y = value, group = 1)) +
  geom_line(color = "#1597A5", size = 1.5) +
  geom_point(color = "#1597A5", size = 4) +
  labs(x = expression(gamma), y = "Relative RMSE") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 25, color = "black", margin = margin(t = 15)), 
    axis.title.y = element_text(size = 25,color = "black"),
    axis.text.x = element_text(size = 20,color = "black"),  
    axis.text.y = element_text(size = 20,color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5))+
    scale_y_continuous(breaks = seq(0.080, max(gama_df$value, na.rm = TRUE), by = 0.01))

ggsave(filename="FigureS6a.png", plot=p1, device="png", dpi=600, units="in", width=8, height=5)


#FigureS6b
df_data <- read.table("~/eNODEconstr/realdata/Synthetic/data/Synthetic_subject_species.txt",header = T,sep="\t")
df_microbe <- df_data[, 4:28]
df_microbe$Sample <- rownames(df_microbe)
df_long <- melt(df_microbe, id.vars = "Sample", variable.name = "Microbe", value.name = "Abundance")
plot_microbe <- ggplot(df_long, aes(x = Microbe, y = Abundance)) +
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
ggsave(filename="FigureS6b.png",plot=plot_microbe,device="png",dpi=600,units="in",width=20,height=10)


#FigureS6c
N <- 25
M <- 4
base_path1 <- "~/eNODEcostr/realdata/Synthetic/results/score/"
true_path <- "/~/eNODEcostr/realdata/Synthetic/data"
setwd(true_path)
abundance_data <- read.table("Synthetic_subject_species.txt",header = T,sep="\t")
microbe_name <- as.data.frame(colnames(abundance_data)[4:(N+3)])
colnames(microbe_name) <- "Microbe_name"
microbe_name$Regulate <- paste0("Microbe",rep(1:N))


setwd(base_path1) 
Score <- read.table("score_normalized_pre_subject_l2.txt",header = T,sep="\t")
Score <- Score %>%
  left_join(microbe_name, by = c("Regulate" = "Regulate")) %>%
  mutate(Regulate = Microbe_name) %>%
  select(-Microbe_name)

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

p1 <-ggplot(result_df, aes(x = "", y = SpearmanCorrelation)) +
  geom_boxplot(fill = "#87CEEB", color = "darkblue", outlier.shape = NA, width = 0.3) +  
  geom_jitter(aes(color = SpearmanCorrelation), width = 0.2, size = 3, shape = 16, alpha = 0.7) + 
  scale_color_gradient(low = "blue", high = "red") +  
  ylab("Spearman Correlation") +
  theme_minimal(base_size = 15) +  
  theme(axis.text.x = element_blank(),  
        axis.title.x = element_blank(),  
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25,  color = "black"),
        legend.position = "none",  
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = "grey80", size = 0.5))
ggsave(filename="FigureS6c.png",plot=p1,device="png",dpi=600,units="in",width=8,height=5)


