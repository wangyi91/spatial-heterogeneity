library(factoextra) # to vis PCA
library(metabaR)
library(cowplot)
library(dplyr)
library(RColorBrewer) # to visualise color palette
library(ggplot2)
library(ggthemes)

source("./00_R-repository/aggregate_to_site.R")

get_pca_plot <- function(dt, color_ind, title, legendpos='none'){
  # function to make PCA plot of a dataset "dt", colored by chosen indicator "color_ind"
  
  color_key = list(long="BurgYl",lat="Emrld", water_depth="PuBuGn")
  
  dt_pca <- prcomp(dt$reads, scale = F)
  
  p <- fviz_pca_ind(dt_pca,
                    col.ind = dt$samples[,color_ind],
                    gradient.cols = hcl.colors(9, color_key[color_ind], rev=T)[c(2:9)], title="",
                    ggtheme=theme_few()) + theme(legend.position = legendpos,
                                                 legend.text = element_text(size=10)) +
    labs(col="water depth\n(m)") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          text = element_text(size=rel(3.5)),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    scale_x_reverse()
  
  p$layers<-p$layers[1] # remove vline and hling at 0s.
  
  return(p)
}


# import data

for (dt_name in c("plants_r2_embl_98","euk_r2_silva_97","cyano_r2_silva_97","cop_r5_EmblSilvaCustom_85")) {
  dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_norm_lgchord", sep=""))
  pa <- dt 
  pa$reads <- (dt$reads>0)*1
  
  
  if (startsWith(dt_name,"p")) {
    a1 <- get_pca_plot(pa,"water_depth","Plant")
  } else if (startsWith(dt_name,"e")) {
    a2 <- get_pca_plot(pa,"water_depth","Eukatyote")
  } else if (startsWith(dt_name,"cy")) {
    a3 <- get_pca_plot(pa,"water_depth","Cyanobacteria")
  } else {
    a4 <- get_pca_plot(pa,"water_depth","Copepod")
  }
  
  
}



# save figure
pdf("ARFig.1_pca_waterdepth.pdf", width=8, height=6)
plot_grid(a1, a2, a3, a4, labels = NULL, nrow=2)
dev.off()

png("ARFig.1_pca_waterdepth.png", width=8, height=6, units = 'in', res = 600)
plot_grid(a1, a2, a3, a4, labels = NULL, nrow=2)
dev.off()








