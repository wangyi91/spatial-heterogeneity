
# Calculate beta diversity decomposed into replacement and nestedness among PCRs from extract in the same site.
# Diversity indices from: Baselga 2010 geb


library(metabaR)
library(dplyr)
library(reshape2)
library(cowplot)
library(ggthemes)
library(ggplot2)

# source this funciton to plot coordinates of samples
source("./05_diversity_analysis/Baselga_2010_geb_appendix_s1.r")


#import data, ASVs appearing in > 1 replicates
dt_name ="plants_r2_embl_98"; colors=c("#F8FCE3", "#B5BE7E");
dt_name ="euk_r2_silva_97"; colors=c("#F8FCE3", "#B5BE7E");
dt_name ="cyano_r2_silva_97"; colors=c("#F0F8FF","#B0C4DE");#colors=c("#F8FCE3", "#B5BE7E");
dt_name ="cop_r5_EmblSilvaCustom_85"; colors=c("#F0F8FF","#B0C4DE");


# import data
dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_2rep", sep=""))
  
if (!startsWith(dt_name, "c")) {
  # aggregate ASVs into taxa (at least to family level) for euk and plant
  dt <- aggregate_motus(dt, groups=dt$motus$TAXID, FUN_agg_motus_sum)
  dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
}
  
dt <- subset_metabarlist(dt, "reads", rowSums(dt$reads)>0)



# presence-absence data
pa <- dt
pa$reads <- (dt$reads>0)*1


# initiate table
beta <- data.frame(0,0,0)[FALSE,]
# create beta-diversity for all sites
for (st in sort(unique(pa$samples$site))) {
  sub <- subset_metabarlist(pa, "samples", pa$samples$site==st)
  beta[st,] <- c(beta.NES(sub$reads),
                 beta.SIM(sub$reads),
                 beta.SOR(sub$reads))
}

names(beta) <- c("beta.NES","beta.SIM","beta.SOR")


beta_long <- 
  reshape2::melt(beta %>% mutate(site=rownames(.)), id.vars="site") %>% 
  filter(variable != "beta.SOR")

dict <- c("plants_r2_embl_98"="Plant", "euk_r2_silva_97"="Eukaryote", "cyano_r2_silva_97"="Cyanobacteria", "cop_r5_EmblSilvaCustom_85"="Copepod")

pdf(paste("Supp.Fig.3c_withinSite_diversity_", substr(dt_name,1,3) ,"_update.pdf", sep=""), width = 7.2, height = 3.6)

beta_long %>%
  ggplot(aes(x=factor(site, levels = rev(rownames(beta))), y=value, col=variable, 
             group=factor(variable))) +
  geom_col(colour="black", size=0.2, width=0.75, aes(fill=variable)) + 
  scale_y_continuous(limits = c(0, 1))+
  xlab("Site") + 
  ylab("Sørensen within-site dissimilarity") + 
  ggtitle(dict[dt_name]) +
  labs(fill='Indices') +
  scale_fill_manual(labels=c("beta.NES\n(nestedness)\n", 
                             "beta.SIM\n(replacement)\n"), values=colors) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  theme(text = element_text(size=14))

dev.off()



