library(tidyverse)
library(ggplot2)

#Figure5

N <- 25
M <- 4

setwd("~/eNODEconstr/realdata/Synthetic/results/trajectory")

#Figure5a
all_microbe_rrmse_mean <- readRDS("all_microbe_rrmse_mean.rds")
all_metabolite_rrmse_mean <- readRDS("all_metabolite_rrmse_mean.rds")
all_microbe_rrmse_mean$Group <- "Microbe"
all_metabolite_rrmse_mean$Group <- "Metabolite"
all_rrmse_mean <- rbind(all_microbe_rrmse_mean,all_metabolite_rrmse_mean)
all_rrmse_mean$Group <- factor(all_rrmse_mean$Group, levels = c("Microbe", "Metabolite"))
all_rrmse_mean$Method <- factor(all_rrmse_mean$Method, levels = c("NODE", "eNODE","eNODEconstr"))

plot_rrmse_mean <- ggplot(all_rrmse_mean, aes(x = Method, y = mean_rrmse, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Group ~ ., scales = "free_y") +  
  scale_fill_manual(values = c("#3F918B","#896191","#df7676")) +
  labs(x = "Method", y = "Relative RMSE", fill = "Method") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 30, color = "black"),
    axis.title.x = element_text(size = 30, color = "black", margin = margin(t = 18)), 
    axis.title.y = element_text(size = 30, color = "black"),
    axis.text.x = element_text(size = 26, color = "black"),  
    axis.text.y = element_text(size = 25.5, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    legend.position = "None",  
    legend.title = element_text(size = 34),  
    legend.text = element_text(size = 34)  
  )

ggsave(filename="Figure5a.png",plot=plot_rrmse_mean,device="png",dpi=600,units="in",width=9,height=11)


#Figure5b
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

#Microbe
microbe_data  <- subset(Score, Group == "Microbe")
medians <- microbe_data %>%
  group_by(Regulate) %>%
  summarize(median_value = median(Score)) %>%
  arrange(median_value)
microbe_data$Regulate <- factor(microbe_data$Regulate, levels = medians$Regulate)
plot1 <- ggplot(microbe_data, aes(x = Score, y = Regulate,fill = Regulate)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.005) +
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
ggsave(filename="Figure5b.png",plot=plot1,device="png",dpi=600,units="in",width=12,height=12)


#Figure5c
metabolite_data  <- subset(Score, Group == "Metabolite")
medians <- metabolite_data %>%
  group_by(Regulate) %>%
  summarize(median_value = median(Score)) %>%
  arrange(median_value)
metabolite_data$Regulate <- factor(metabolite_data$Regulate, levels = medians$Regulate)
plot2 <- ggplot(metabolite_data, aes(x = Score, y = Regulate,fill = Regulate)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.01) +
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
ggsave(filename="Figure5c.png",plot=plot2,device="png",dpi=600,units="in",width=12,height=12)


#Figure5d
setwd("~/eNODEcostr/realdata/Synthetic/results/interaction")
mu_pre <- readRDS("mu_pre.rds") 
mu_true_binary <- readRDS("~/eNODEcostr/realdata/Synthetic/data/Input/mu.rds")

mu_pre_only <- mu_pre 
mu_pre_only[!(mu_true_binary == 0 & mu_pre != 0)] <- NA 

mu_pre_only <- as.data.frame(mu_pre_only)
mu_pre_only_long <- mu_pre_only %>%
  rownames_to_column(var = "Species") %>%
  gather(key = "Metabolite", value = "Value", -Species) %>%
  filter(Value != "NA" )  

edges <- mu_pre_only_long %>%
  select(Species, Metabolite, Value)

nodes <- data.frame(name = unique(c(edges$Species, edges$Metabolite)),
                    type = ifelse(unique(c(edges$Species, edges$Metabolite)) %in% edges$Species, "Species", "Metabolite"))

nodes$type <- factor(nodes$type, levels = c("Species", "Metabolite"))

graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
V(graph)$type <- ifelse(V(graph)$name %in% edges$Species, TRUE, FALSE)
edge_count <- ecount(graph)
set.seed(12)  
layout_circle <- create_layout(graph, layout = "circle")
layout_circle$color <- ifelse(layout_circle$type == TRUE, "#1597A5", "#ffbe7a")

plot_circle1 <- ggraph(layout_circle) + 
  geom_edge_link(aes(edge_colour = ifelse(Value > 0, "red", "blue"),
                     edge_width = abs(Value)),
                 edge_alpha = 0.15) + 
  geom_node_point(colour = layout_circle$color, size = 12) +  
  geom_node_text(data = subset(layout_circle, color == "#ffbe7a"), 
                 aes(label = name, 
                     x = x + 0.15,  
                     y = y),       
                 size = 6, 
                 repel = FALSE) +  
  geom_node_text(data = subset(layout_circle, color != "#ffbe7a"), 
                 aes(label = name), 
                 size = 6, 
                 repel = TRUE) +  
  scale_edge_colour_identity() +  
  theme_void() +
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, -0.2), "cm"))
ggsave(filename="Figure5d.png", plot=plot_circle1, device="png", dpi=600, units="in", width=12.5, height=10)


#Figure5e
mu_pre_shared <- mu_pre 
mu_pre_shared[!(mu_true_binary != 0 & mu_pre != 0)] <- NA 
mu_pre_shared <- as.data.frame(mu_pre_shared)
mu_pre_shared_long <- mu_pre_shared %>%
  rownames_to_column(var = "Species") %>%
  gather(key = "Metabolite", value = "Value", -Species) %>%
  filter(Value != "NA" )  

edges <- mu_pre_shared_long %>%
  select(Species, Metabolite, Value)

nodes <- data.frame(name = unique(c(edges$Species, edges$Metabolite)),
                    type = ifelse(unique(c(edges$Species, edges$Metabolite)) %in% edges$Species, "Species", "Metabolite"))
nodes$type <- factor(nodes$type, levels = c("Species", "Metabolite"))
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
V(graph)$type <- ifelse(V(graph)$name %in% edges$Species, TRUE, FALSE)
edge_count <- ecount(graph)
set.seed(12345)  
layout_circle <- create_layout(graph, layout = "circle")
layout_circle$color <- ifelse(layout_circle$type == TRUE, "#1597A5", "#ffbe7a")
plot_circle2 <- ggraph(layout_circle) + 
  geom_edge_link(aes(edge_colour = ifelse(Value > 0, "red", "blue"),
                     edge_width = abs(Value)),
                 edge_alpha = 0.15) +  
  geom_node_point(colour = layout_circle$color, size = 12) +  
  geom_node_text(data = subset(layout_circle, color == "#ffbe7a"), 
                 aes(label = name, 
                     x = x + 0.15,  
                     y = y),       
                 size = 6, 
                 repel = FALSE) +  
  geom_node_text(data = subset(layout_circle, color != "#ffbe7a"), 
                 aes(label = name), 
                 size = 6, 
                 repel = TRUE) +  
  scale_edge_colour_identity() +  
  theme_void() +
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, -0.3), "cm"))
ggsave(filename="Figure5e.png", plot=plot_circle2, device="png", dpi=600, units="in", width=12.5, height=10)
