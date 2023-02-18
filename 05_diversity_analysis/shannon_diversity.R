# this script calculates Shannon alpha diversity based on taxa abundance

library(ade4)
library(vegan)
library(metabaR)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggthemes)
library(RColorBrewer)



for (dt_name in c("plants_r2_embl_98","euk_r2_silva_97","cyano_r2_silva_97","cop_r5_EmblSilvaCustom_85")) {
  
  # load data
  dt <- readRDS(paste("./00_R-repository/",dt_name, "_agg_2rep_motuAgg_norm", sep=""))
  
  if (dt_name == "euk_r2_silva_97"){
    dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  } else if (dt_name == "plants_r2_embl_98"){
    dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  } else if (dt_name == "cyano_r2_silva_97"){
    dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  } else {
    dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  }

  dt <- subset_metabarlist(dt, "reads", rowSums(dt$reads)>0)
  
  
  # make a table with all info needed
  alphadiv <- diversity(dt$reads, index="shannon") %>% exp %>% as.data.frame %>% `colnames<-`("shannon") %>% 
    mutate(sample_id=rownames(.)) %>% 
    left_join(., dt$sample %>% 
              mutate(sample_id=rownames(.)) %>% 
              select(sample_id, sediment, long, lat, water_depth, site), 
              by="sample_id") %>%
    group_by(site) %>% mutate(same_core=as.integer(factor(sediment)))
  
  
  # save to each dataset separately
  if (startsWith(dt_name,"p")){
    a1<-alphadiv
  } else if (startsWith(dt_name,"eu")){
    a2<-alphadiv
  } else if (startsWith(dt_name,"cy")){
    a3<-alphadiv
  } else {
    a4<-alphadiv
  }
}


pdf("Fig.1b_shannon_diversity.pdf", width = 7.2, height = 4.8)
# make plot combining all four datasets
ggplot() +
  geom_point(data=a1, aes(x=site, y=shannon), color=brewer.pal(n = 8, name = "Dark2")[6], alpha=0.8,size=3, shape=16) +
  geom_point(data=a2, aes(x=site, y=shannon), color=brewer.pal(n = 8, name = "Dark2")[3], alpha=0.8,size=3, shape=17) +
  geom_point(data=a3, aes(x=site, y=shannon), color="black", alpha=0.5,size=3.3, shape=18) +
  geom_point(data=a4, aes(x=site, y=shannon), color=brewer.pal(n = 8, name = "Dark2")[4], stroke=0.9,alpha=0.8,size=2, shape=6) +
  scale_x_discrete(limits=rev) +
  scale_y_log10(limits=c(1,100), breaks = c(1, 3, 10,30,100)) +
  annotation_logticks(sides = "l",short = unit(0.1, "cm"),
                      mid = unit(0.15, "cm"),
                      long = unit(0.2, "cm"),size=0.4) +
  theme_classic() + 
  ylab("Shannon diversity of taxa") + 
  xlab("Site") + 
  theme(axis.text.x = element_text(angle = 45,hjust=0.5, vjust=0.5, size=12),
        text = element_text(size=rel(4.5)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)))#,
        #panel.border = element_rect(colour = "black", fill=NA, size=1.4))  
dev.off()




##################################################################################
library(rlang)
library(ggpmisc)
dict <- c("long"="Longitude (°E)", "lat"="Latitude (°N)", "water_depth"="Water depth (m)")

plot_alpha <- function(data, x, y, Color, Shape, lm=FALSE, labelx){
  gg <- ggplot(data, aes({{x}},{{y}})) + 
      geom_point(color=Color, size=3, shape=Shape) +
      xlab(dict[deparse(substitute(x))])
  
    if (lm==TRUE){
      gg <- gg +
        geom_smooth(method='lm', size=0.7, color="black", linetype="dashed") +
        stat_poly_eq(aes(label = paste(after_stat(rr.label), sep="*\", \"*")), label.x = labelx) +
        stat_poly_eq(aes(label = paste("italic(p)-value", formatC(after_stat(p.value), format = "e", digits = 1), 
                                     sep="*\" = \"*")),
                     label.x = labelx, label.y = 0.8) + 
        theme_classic() +
        theme(text = element_text(size=rel(4)),
              axis.title.x = element_text(size=15),
              axis.title.y = element_blank())
      } else {
        gg <- gg + theme_classic() +
          theme(text = element_text(size=rel(4)),
              axis.title.x = element_text(size=15),
              axis.title.y = element_blank())
      } 
  return(gg)
}


long1 <- plot_alpha(a1, long, shannon, "black", 16)
lat1 <- plot_alpha(a1, lat, shannon, "black", 17)
wd1 <- plot_alpha(a1, water_depth, shannon, "black", 18)

long2 <- plot_alpha(a2, long, shannon, "black", 16)
lat2 <- plot_alpha(a2, lat, shannon, "black", 17)
wd2 <- plot_alpha(a2, water_depth, shannon, "black", 18)

long3 <- plot_alpha(a3, long, shannon, "black", 16, lm=TRUE, labelx=0.1)
lat3 <- plot_alpha(a3, lat, shannon, "black", 17, lm=TRUE, labelx=0.3)
wd3 <- plot_alpha(a3, water_depth, shannon, "black", 18)

long4 <- plot_alpha(a4, long, shannon, "black", 16)
lat4 <- plot_alpha(a4, lat, shannon, "black", 17)
wd4 <- plot_alpha(a4, water_depth, shannon, "black", 18)

pdf("Supp.Fig.3a_shannon_diversity.pdf", width = 12, height = 6.2)
plot_grid(long1, long2, long3, long4,
          lat1, lat2, lat3, lat4,
          wd1, wd2, wd3, wd4, labels = NULL, nrow=3)
dev.off()




