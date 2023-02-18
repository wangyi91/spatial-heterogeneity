# This script calculates beta diversity (Sørensen dissimilarity) between two DNA extracts using presence/absence data
# Diversity indices from: Legendre 2014

library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(metabaR)
source("./05_diversity_analysis/beta.div.comp_Legendre_2014_geb.R")



#### Part I: Make Fig.1c ####

for (dt_name in c("plants_r2_embl_98","euk_r2_silva_97")) { #,"cyano_r2_silva_97","cop_r5_EmblSilvaCustom_85"
  # import data
  dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_2rep", sep=""))
  
  # aggregate ASVs into taxa
  dt <- aggregate_motus(dt, groups=dt$motus$TAXID, FUN_agg_motus_sum)
  
  if (dt_name == "euk_r2_silva_97"){
    dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  } else if (dt_name == "plants_r2_embl_98"){
    dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  }
  dt <- subset_metabarlist(dt, "reads", rowSums(dt$reads)>0)
  
  
  # presence-absence data
  pa <- dt
  pa$reads <- (dt$reads>0)*1
  
  # calculate betas for whole data
  mat_beta <- beta.div.comp(pa$reads,coef="BS", quant=F)
  
  # dist to data.frame for pairwise data
  repl <- reshape2::melt(as.matrix(mat_beta$repl), varnames = c("row", "col"), value.name = "repl") %>% filter(repl>0)
  
  rich <- reshape2::melt(as.matrix(mat_beta$rich), varnames = c("row", "col"), value.name = "rich") %>% filter(rich>0)
  
  D <- reshape2::melt(as.matrix(mat_beta$D), varnames = c("row", "col"), value.name = "D") %>% filter(D>0)
  
  betas <- left_join(repl, rich, by=c("row","col")) %>% left_join(., D, by=c("row","col"))
  
  # calculate long/lat distances between two samples
  betas$dist.long <- abs(pa$samples[betas$row,]$long - pa$samples[betas$col,]$long) * 111.321
  betas$dist.lat <- abs(pa$samples[betas$row,]$lat - pa$samples[betas$col,]$lat) * 111.321
  betas$dist.geo <- sqrt(betas$dist.long^2 + betas$dist.lat^2)
  
  # save to each dataset separately
  if (startsWith(dt_name,"p")){
    b1<-betas
  } else if (startsWith(dt_name,"eu")){
    b2<-betas
  } else if (startsWith(dt_name,"cy")){
    b3<-betas
  } else {
    b4<-betas
  }
}


pdf("Fig.1c_beta_diversity.pdf", width = 7.2, height = 5)

ggplot() + 
  #geom_point(data=b3, aes(x=dist.geo, y=D), color="black",size=1.2, shape=5, alpha=0.2) +
  #geom_point(data=b4, aes(x=dist.geo, y=D), color=brewer.pal(n = 8, name = "Dark2")[4],size=1.2, shape=25, alpha=0.16) +
  geom_point(data=b2, aes(x=dist.geo, y=D), color=brewer.pal(n = 8, name = "Dark2")[3],size=1.6, shape=17, alpha=0.1) +
  geom_smooth(data=b2,aes(x=dist.geo, y=D), method='lm', size=0.7, color="black", linetype="dotdash") +
  geom_point(data=b1, aes(x=dist.geo, y=D), color=brewer.pal(n = 8, name = "Dark2")[6],size=1.6, shape=16, alpha=0.1) + 
  geom_smooth(data=b1, aes(x=dist.geo, y=D), method='lm', size=0.7, color="black", linetype="dashed") +
  xlab("Geographic distance (km)") + 
  ylab("Across-site Sørensen dissimilarity") +
  ylim(0, 1) + 
  theme_classic() +
  theme(text = element_text(size=rel(4.5)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)))#,
        #panel.border = element_rect(colour = "black", fill=NA, size=1.4))
dev.off()

fit1 = lm(D~dist.geo, data=b1); summary(fit1)
fit2 = lm(D~dist.geo, data=b2); summary(fit2)



#### Part II: Make Supp.Fig.3b ####

library(rlang)
library(ggpmisc)

plot_beta <- function(data, x, y, Color, Shape){
  ggplot(data, aes({{x}},{{y}})) + 
    geom_point(color=Color, size=1.5, shape=Shape, alpha=0.1) +
    geom_smooth(method='lm', size=0.7, color="black", linetype="dashed") +
    #xlab("Geographic distance (km)") + 
    #ylab("Within-site Sørensen dissimilarity") +
    ylim(0, 1) + 
    theme_classic() +
    theme(text = element_text(size=rel(4.5)),
          #axis.title.y = element_text(margin = margin(r = 10)),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                   after_stat(rr.label), 
                                   sep="*\", \"*"))) +
    stat_poly_eq(aes(label = paste("italic(p)-value", 
                                   formatC(after_stat(p.value), format = "e", digits = 1), sep="*\" = \"*")),
                     label.y = 0.85)
}

for (dt_name in c("plants_r2_embl_98","euk_r2_silva_97","cyano_r2_silva_97","cop_r5_EmblSilvaCustom_85")) { 
  # import data
  dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_2rep", sep=""))
  
  if (!startsWith(dt_name, "c")){
    # aggregate ASVs into taxa for euk and plant
    dt <- aggregate_motus(dt, groups=dt$motus$TAXID, FUN_agg_motus_sum)
  }
  
  if (dt_name == "euk_r2_silva_97"){
    dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  } else if (dt_name == "plants_r2_embl_98"){
    dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  } #else if (dt_name == "cyano_r2_silva_97"){
    #dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
  #}
  
  dt <- subset_metabarlist(dt, "reads", rowSums(dt$reads)>0)
  
  # presence-absence data
  pa <- dt
  pa$reads <- (dt$reads>0)*1
  
  # calculate betas for whole data
  mat_beta <- beta.div.comp(pa$reads,coef="BS", quant=F)
  
  # dist to data.frame for pairwise data
  repl <- reshape2::melt(as.matrix(mat_beta$repl), varnames = c("row", "col"), value.name = "repl") %>% filter(repl>0)
  
  rich <- reshape2::melt(as.matrix(mat_beta$rich), varnames = c("row", "col"), value.name = "rich") %>% filter(rich>0)
  
  D <- reshape2::melt(as.matrix(mat_beta$D), varnames = c("row", "col"), value.name = "D") %>% filter(D>0)
  
  betas <- left_join(repl, rich, by=c("row","col")) %>% left_join(., D, by=c("row","col"))
  
  # calculate long/lat distances between two samples
  betas$dist.long <- abs(pa$samples[betas$row,]$long - pa$samples[betas$col,]$long) * 111.321
  betas$dist.lat <- abs(pa$samples[betas$row,]$lat - pa$samples[betas$col,]$lat) * 111.321
  betas$dist.geo <- sqrt(betas$dist.long^2 + betas$dist.lat^2)
  
  # save to each dataset separately
  if (startsWith(dt_name,"p")){
    D1 <- plot_beta(betas, dist.geo, D, "black", 16)
    repl1 <- plot_beta(betas, dist.geo, repl, "black", 16)
    rich1 <- plot_beta(betas, dist.geo, rich, "black", 16)
  } else if (startsWith(dt_name,"eu")){
    D2 <- plot_beta(betas, dist.geo, D, "black", 16)
    repl2 <- plot_beta(betas, dist.geo, repl, "black", 16)
    rich2 <- plot_beta(betas, dist.geo, rich, "black", 16)
  } else if (startsWith(dt_name,"cy")){
    D3 <- plot_beta(betas, dist.geo, D, "royal blue", 17)
    repl3 <- plot_beta(betas, dist.geo, repl,"royal blue", 17)
    rich3 <- plot_beta(betas, dist.geo, rich,"royal blue", 17)
  } else {
    D4 <- plot_beta(betas, dist.geo, D, "royal blue", 17)
    repl4 <- plot_beta(betas, dist.geo, repl, "royal blue", 17)
    rich4 <- plot_beta(betas, dist.geo, rich, "royal blue", 17)
  }
}

pdf("Supp.Fig.3b_beta_diversity.pdf", width = 13.5, height = 9)

plot_grid(D1, D2, D3, D4,
          repl1, repl2, repl3, repl4,
          rich1, rich2, rich3, rich4, labels = NULL, nrow=3)

dev.off()

