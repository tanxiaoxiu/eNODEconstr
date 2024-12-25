library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(igraph)
library(reticulate)
library(tidyr)


############Figure6
base_path1 <- "~/eNODEconstr/realdata/Fiber/results/trajectory"
true_path <- "~/eNODEconstr/realdata/Fiber/data"
setwd(true_path)
abundance_data <- read.table("filter_fiber_SCFA.txt",header = T,sep="\t")
microbe_name <- as.data.frame(colnames(abundance_data)[6:(N+5)])
colnames(microbe_name) <- "Microbe_name"
microbe_name$Regulate <- paste0("Microbe",rep(1:N))

#Figure6a
path <- paste0(base_path, "/Rs")
setwd(path)
#####Prediction
Score <- read.table("score_normalized_pre_subject_l2.txt",header = T,sep="\t")
Score <- Score %>%
  left_join(microbe_name, by = c("Regulate" = "Regulate")) %>%
  mutate(Regulate = Microbe_name) %>%
  select(-Microbe_name)

#Microbe
microbe_data  <- subset(Score, Group == "Microbe")
microbe_data <- na.omit(microbe_data)
medians <- microbe_data %>%
  group_by(Regulate) %>%
  summarize(median_value = median(Score)) %>%
  arrange(median_value)
microbe_data$Regulate <- factor(microbe_data$Regulate, levels = medians$Regulate)

plot1 <- ggplot(microbe_data, aes(x = Score, y = Regulate,fill = Regulate)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.02) +
  scale_fill_manual(values = rep("#1597A5", length(unique(microbe_data$Regulate)))) +
  labs(x = NULL, y = "Microbe") +
  xlim(0, NA) +
  theme_ridges() +
  theme(axis.title.y = element_blank(), 
        legend.position = "none",
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey90"),
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25))
plot1
ggsave(filename="Figure6a.png",plot=plot1,device="png",dpi=600,units="in",width=12,height=12)

#Figure6b
metabolite_data  <- subset(Score, Group == "Metabolite")
metabolite_data <- na.omit(metabolite_data)
medians <- metabolite_data %>%
  group_by(Regulate) %>%
  summarize(median_value = median(Score)) %>%
  arrange(median_value)
metabolite_data$Regulate <- factor(metabolite_data$Regulate, levels = medians$Regulate)

plot2 <- ggplot(metabolite_data, aes(x = Score, y = Regulate,fill = Regulate)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.02) +
  scale_fill_manual(values = rep("#ffbe7a", length(unique(metabolite_data$Regulate)))) +
  labs(x = NULL, y = "Microbe") +
  xlim(0, NA) +
  theme_ridges() +
  theme(axis.title.y = element_blank(), 
        legend.position = "none",
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey90"),
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25))
plot2
ggsave(filename="Figure6b.png",plot=plot2,device="png",dpi=600,units="in",width=12,height=12)

#Figure6c
path <- paste0(base_path, "/In")
setwd(path)
Score <- read.table("score_normalized_pre_subject_l2.txt",header = T,sep="\t")
Score <- Score %>%
  left_join(microbe_name, by = c("Regulate" = "Regulate")) %>%
  mutate(Regulate = Microbe_name) %>%
  select(-Microbe_name)

#Microbe
microbe_data  <- subset(Score, Group == "Microbe")
microbe_data <- na.omit(microbe_data)
medians <- microbe_data %>%
  group_by(Regulate) %>%
  summarize(median_value = median(Score)) %>%
  arrange(median_value)
microbe_data$Regulate <- factor(microbe_data$Regulate, levels = medians$Regulate)

plot1 <- ggplot(microbe_data, aes(x = Score, y = Regulate,fill = Regulate)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.008) +
  scale_fill_manual(values = rep("#1597A5", length(unique(microbe_data$Regulate)))) +
  labs(x = NULL, y = "Microbe") +
  xlim(0, NA) +
  theme_ridges() +
  theme(axis.title.y = element_blank(), 
        legend.position = "none",
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey90"),
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25),
        plot.margin = margin(t = 5.5, r = 17.5, b = 5.5, l = 5.5, unit = "pt"))
ggsave(filename="Figure6c.png",plot=plot1,device="png",dpi=600,units="in",width=12,height=12 )

#Metabolite
metabolite_data  <- subset(Score, Group == "Metabolite")
metabolite_data <- na.omit(metabolite_data)
medians <- metabolite_data %>%
  group_by(Regulate) %>%
  summarize(median_value = median(Score)) %>%
  arrange(median_value)
metabolite_data$Regulate <- factor(metabolite_data$Regulate, levels = medians$Regulate)

plot2 <- ggplot(metabolite_data, aes(x = Score, y = Regulate,fill = Regulate)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.008) +
  scale_fill_manual(values = rep("#ffbe7a", length(unique(metabolite_data$Regulate)))) +
  labs(x = NULL, y = "Microbe") +
  xlim(0, NA) +
  theme_ridges() +
  theme(axis.title.y = element_blank(), 
        legend.position = "none",
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey90"),
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25),
        plot.margin = margin(t = 5.5, r = 17.5, b = 5.5, l = 5.5, unit = "pt"))
ggsave(filename="Figure6d.png",plot=plot2,device="png",dpi=600,units="in",width=12,height=12)




