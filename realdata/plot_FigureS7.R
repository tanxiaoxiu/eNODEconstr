library(dplyr)
library(reshape2)
library(ggplot2)

###########FigureS7
#FigureS7a
setwd("~/eNODEconstr/realdata/Fiber/results/trajectory")
gama <- read.table("gama.txt",header = T,sep="\t")
colnames(gama) <- c("Group", "1e0","1e-1", "1e-2","1e-3","1e-4","1e-5","1e-6")
gama$Group <- factor(gama$Group, levels = c("Con", "Rs", "In"))
gama_long <- gather(gama, key = "gama", value = "value", -Group)
gama_long$gama <- factor(gama_long$gama, levels = c("1e0","1e-1", "1e-2","1e-3","1e-4","1e-5","1e-6"))
custom_colors <- c("Con" = "#1597A5", "Rs" = "#FEB3AE", "In" = "#FFC24B")
p1 <- ggplot(gama_long, aes(x = gama, y = value, group = Group, color = Group)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  scale_color_manual(values = custom_colors) +
  labs(x = expression(gamma), y = "Relative RMSE") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 25, color = "black", margin = margin(t = 15)), 
    axis.title.y = element_text(size = 25,color = "black"),
    axis.text.x = element_text(size = 20,color = "black"),  
    axis.text.y = element_text(size = 20,color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 22))

ggsave(filename="FigureS7a.png",plot=p1,device="png",dpi=600,units="in",width=10,height=5)


#FigureS7b
###########mean
setwd("~/eNODEconstr/realdata/Fiber/results/trajectory")
all_rrmse_mean <- readRDS("all_rrmse_mean.rds")
all_rrmse_mean$Method <- factor(all_rrmse_mean$Method, levels = c("NODE", "eNODE","eNODEconstr"))

all_rrmse_mean$mean_rrmse <- as.numeric(as.character(all_rrmse_mean$mean_rrmse))
all_rrmse_mean$Group <- factor(all_rrmse_mean$Group, levels = c("Control", "Rs", "In"))
all_rrmse_mean$Method <- factor(all_rrmse_mean$Method, levels = c("NODE", "eNODE","eNODEconstr"))

g <- ggplot(all_rrmse_mean, aes(x=Group, y=mean_rrmse, fill=Method))
p2 <- g + geom_col(position = position_dodge()) +
  scale_fill_manual(values = c("#3F918B","#896191","#df7676")) +
  xlab("Group") +
  ylab("Relative RMSE") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +  
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30)) +
  theme(panel.grid.major = element_line(colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.background = element_rect(fill="transparent", colour=NA),
        panel.grid.minor = element_blank()) +
  theme(legend.spacing.y = unit(0.5, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE))

ggsave(filename="FigureS7b.png",plot=p2,device="png",dpi=600,units="in",width=14,height=8)

