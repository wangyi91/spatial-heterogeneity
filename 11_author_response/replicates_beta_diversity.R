# This script calculates beta diversity (Sørensen dissimilarity) between two DNA extracts 
# using taxon ASV count and taxon presence/absence data for the eukaryote dataset

# Distances between samples collected from separate cores at the same site were estimated to be 20 cm apart 
# samples from the surface sediment of the same core were estimated to be 2 cm apart.

# Diversity indices from: Legendre 2014

library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpmisc)
library(metabaR)
source("./05_diversity_analysis/beta.div.comp_Legendre_2014_geb.R")


# import data
dt <- readRDS(paste("./00_R-repository/euk_r2_silva_97_agg_2rep_motuAgg_norm", sep=""))
  
# keep taxa identified at least to family
dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$family_name))
dt <- subset_metabarlist(dt, "reads", rowSums(dt$reads)>0)

# choose using ASV count or presence/absence data
tr <- dt
trans="presence/absence"; fig_no="a"; Color="cornflowerblue"
trans="ASV count"; fig_no="b"; Color="darkmagenta"

if (trans=="presence/absence"){
  tr$reads <- (dt$reads>0)*1
  # calculate betas for whole data
  mat_beta <- beta.div.comp(tr$reads,coef="BS", quant=F)
} else {
  # calculate betas for whole data
  mat_beta <- beta.div.comp(tr$reads,coef="BS", quant=T)
  
}

  
  
# dist to data.frame for pairwise data
repl <- reshape2::melt(as.matrix(mat_beta$repl), varnames = c("row", "col"), value.name = "repl") %>% filter(repl>0)
  
rich <- reshape2::melt(as.matrix(mat_beta$rich), varnames = c("row", "col"), value.name = "rich") %>% filter(rich>0)
  
D <- reshape2::melt(as.matrix(mat_beta$D), varnames = c("row", "col"), value.name = "D") %>% filter(D>0)
  
betas <- left_join(repl, rich, by=c("row","col")) %>% left_join(., D, by=c("row","col"))
  
# calculate long/lat distances between two samples
betas$dist.long <- abs(tr$samples[betas$row,]$long - tr$samples[betas$col,]$long) * 111.321
betas$dist.lat <- abs(tr$samples[betas$row,]$lat - tr$samples[betas$col,]$lat) * 111.321
betas$dist.geo <- sqrt(betas$dist.long^2 + betas$dist.lat^2)

# subset data
bsub <- betas %>% select(-dist.long, -dist.lat) %>% filter(dist.geo < 1) %>% 
  mutate(dist.geo.new = ifelse(tr$samples[.$row,]$sediment==tr$samples[.$col,]$sediment, 2e-5,
                           ifelse(tr$samples[.$row,]$site==tr$samples[.$col,]$site, 2e-4,
                                  dist.geo)))


pdf(paste("Author_response_figure_1", fig_no, ".pdf", sep=""), width = 7.2, height = 5)

ggplot(data=bsub, aes(x=dist.geo.new, y=D)) + 
  geom_point(color=Color, size=1.6, shape=16, alpha=0.2) +
  geom_smooth(method='lm', size=0.7, color="black", linetype="dotdash") +
  xlab("Geographic distance (km)") + 
  ylab(paste("Across-site Sørensen dissimilarity\ntaxon ", trans, sep="")) +
  ylim(0, 1) + 
  theme_classic() +
  theme(text = element_text(size=rel(4.5)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                               after_stat(rr.label), 
                               sep="*\", \"*"))) +
  stat_poly_eq(aes(label = paste("italic(p)-value", 
                                 formatC(after_stat(p.value), format = "e", digits = 1), sep="*\" = \"*")),
               label.y = 0.85)

dev.off()


