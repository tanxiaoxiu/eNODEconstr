library(tidyverse)
library(ggplot2)
library(patchwork)


############plot FigureS5
N <- 20
M <- 20
compare_path <- "~/eNODEconstr/simulation/comparison/DKI"

path <- paste0(compare_path, "/n", N, "m", M)
setwd(path)
all_results <- readRDS("all_results_DKI.rds")

all_score_r_rmse_mean <- all_results$all_score_r_rmse_mean
all_score_r_rmse_mean <- na.omit(all_score_r_rmse_mean)

g <- ggplot(all_score_r_rmse_mean, aes(x=as.character(SubjectN)))
plot1.1 <- g + geom_col(aes(fill = Method,y= r_rmse), position = position_dodge())+
  scale_fill_manual(values = c("#82b0d2","#df7676"))+
  xlab(expression(italic(S)))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
p1.1 <-  plot1.1 + theme(legend.spacing.y = unit(1, 'cm')) + guides(fill = guide_legend(byrow = TRUE))

all_score_spearman_mean <- all_results$all_score_spearman_mean
all_score_spearman_mean <- na.omit(all_score_spearman_mean)
g <- ggplot(all_score_spearman_mean, aes(x=as.character(SubjectN)))
plot1.2 <- g + geom_col(aes(fill = Method,y= spearman), position = position_dodge())+
  scale_fill_manual(values = c("#82b0d2","#df7676"))+
  xlab(expression(italic(S)))+
  ylab("Spearman correlation")+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
p1.2 <-  plot1.2 + theme(legend.spacing.y = unit(1, 'cm')) + guides(fill = guide_legend(byrow = TRUE))

plot_all <- p1.1 + p1.2 + 
  plot_layout(ncol = 2, nrow = 1) + 
  plot_annotation(tag_levels = 'a')
plot_all <- plot_all + plot_layout(guides = 'collect')
ggsave(filename="FigureS5.png",plot=plot_all,device="png",dpi=600,units="in",width=20,height=8)




