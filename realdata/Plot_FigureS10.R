library(pheatmap)
library(gridExtra)
library(grid)

######FigureS10
stewd("~/eNODEconstr/realdata/Fiber/results/interaction/Rs") 
mu_pre <- readRDS("mu_pre.rds")
p1 <- pheatmap(mu_pre,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "none",
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-1.5, 1.5, length.out = 101),  
               center = 0,
               fontsize_row = 13,        
               fontsize_col = 13,        
               legend = TRUE,
               main = "")

ggsave(filename="FigureS10a.png",plot=p1,device="png",dpi=600,units="in",width=5.5,height=6)

stewd("~/eNODEconstr/realdata/Fiber/results/interaction/In") 
mu_pre <- readRDS("mu_pre.rds")
p2 <- pheatmap(mu_pre,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "none",
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-1.5, 1.5, length.out = 101),  
               center = 0,
               fontsize_row = 13,        
               fontsize_col = 13,        
               legend = TRUE,
               main = "")

ggsave(filename="FigureS10b.png",plot=p2,device="png",dpi=600,units="in",width=5.5,height=6)
