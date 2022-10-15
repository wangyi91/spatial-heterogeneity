# this script calculatea Shannon alpha diversity based on taxa abundance, then generates Fig.1C

library(ade4)
library(vegan)
library(metabaR)
library(cowplot)
library(dplyr)
library(ggthemes)
library(RColorBrewer)



for (dt_name in c("plants_r2_embl_98","euk_r2_silva_97","cyano_r2_silva_97","cop_r5_EmblSilvaCustom_85")) {
  
  # load data
  dt<-readRDS(paste("./00_R-repository/",dt_name, "_agg_2rep_motuAgg_norm", sep=""))
  
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


pdf("Fig.1C_shannon_diversity.pdf", width = 7.2, height = 4.8)
# make plot combining all four datasets
ggplot() +
  geom_point(data=a1, aes(x=site, y=shannon), color=brewer.pal(n = 8, name = "Dark2")[6], alpha=0.8,size=3, shape=16) +
  geom_point(data=a2, aes(x=site, y=shannon), color=brewer.pal(n = 8, name = "Dark2")[3], alpha=0.8,size=3, shape=17) +
  geom_point(data=a3, aes(x=site, y=shannon), color="black", alpha=0.5,size=3.3, shape=18) +
  geom_point(data=a4, aes(x=site, y=shannon), color=brewer.pal(n = 8, name = "Dark2")[4], stroke=0.9,alpha=0.8,size=2, shape=6) +
  scale_x_discrete(limits=rev) +
  scale_y_log10(limits=c(1,100), breaks = c(1, 3, 10,30,100))+
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
