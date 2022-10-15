library(FactoMineR)
library(factoextra) # to vis PCA
library(metabaR)
library(cowplot)
library(dplyr)


library(RColorBrewer) # to visualise color palette

library(ggplot2)
library(ggthemes)


get_pca_plot <- function(dt, color_ind, title){
  # function to make PCA plot of a dataset "dt", colored by chosen indicator "color_ind"
  
  color_key = list(long="BurgYl",lat="Emrld", water_depth="Oslo")
  
  dt_pca <- prcomp(dt$reads, scale = F)
  
  p <- fviz_pca_ind(dt_pca,
                    col.ind = dt$samples[,color_ind],
                    label = 'none',
                    gradient.cols = hcl.colors(9, color_key[color_ind], rev=F), title="",
                    ggtheme=theme_few()) + theme(legend.position = "none") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          text = element_text(size=rel(4)),
          #axis.title.y = element_text(margin = margin(r = 10)),
          #axis.title.x = element_text(margin = margin(t = 10)),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_reverse()
  
  p$layers<-p$layers[1] # remove vline and hling at 0s.
  
  return(p)
}


#import data, already motu-aggregated and log-chord transformed
dt1 <- readRDS("./00_R-repository/plants_r2_embl_98_agg_norm_lgchord")

dt2 <- readRDS("./00_R-repository/euk_r2_silva_97_agg_norm_lgchord")

dt3 <- readRDS("./00_R-repository/cyano_r2_silva_97_agg_norm_lgchord")

dt4 <- readRDS("./00_R-repository/cop_r2_EmblSilvaCustom_85_agg_norm_lgchord")





a1 <- get_pca_plot(dt1,"long","PLANT")
a2 <- get_pca_plot(dt2,"long","EUK")
a3 <- get_pca_plot(dt3,"long","CYA")
a4 <- get_pca_plot(dt4,"long","COP")




# save longitude plot without legend
pdf("Fig.1B_pca_long_noLegend.pdf", width=12, height=2.55)

plot_grid(a1, a2, a3, a4, labels = NULL, nrow=1)

dev.off()








