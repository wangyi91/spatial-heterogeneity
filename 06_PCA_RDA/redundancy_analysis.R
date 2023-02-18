# RDA using longitude, latitude and sedimentological data as explanatory variables

library(vegan) # function rda
library(metabaR)
library(dplyr)
library(fuzzyjoin)
library(cowplot)
library(shape) # to use the Arrows function

# import data. Choose one
dt_name = "plants_r2_embl_98"
dt_name = "euk_r2_silva_97"
dt_name = "cyano_r2_silva_97"
dt_name = "cop_r5_EmblSilvaCustom_85"

dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_norm_lgchord", sep=""))



# import sedimentological data table
sedi_dt <- read.table("./06_PCA_RDA/surfsedi_sedimentological_data.txt", header=T, sep='\t')

# match sedi data to 'samples' table, allowing 100 m max distance mismatch
sedi <- difference_left_join(dt$samples, sedi_dt, by=c("UTM.32T.E","UTM.32T.N"),  max_dist = 100, distance_col="dist") %>%
  mutate(UTM.32T.E=UTM.32T.E.x, UTM.32T.N=UTM.32T.N.x)


# remove non-numeric and useless columns
sedi[,c("sediment","KERN","VCoarseSand","CoarseSand","GoodnessofFit","Kies",
        "UTM.32T.N.x","UTM.32T.N.y","UTM.32T.E.x","UTM.32T.E.y","dist")] <- NULL


# scale numeric variables
sedi_st <- scale(sedi %>% select(-c(site)), center = T, scale = T) %>% as.data.frame



# Model 1: long, lat and water depth as explanatory variables.
mod1 <- rda(dt$reads ~ long + lat + water_depth, sedi_st, na.action=na.exclude)
RsquareAdj(mod1)$adj.r.squared

summary(mod1)
#RsquareAdj(mod1)$r.squared

# permutation test of explanatory factors
anova.cca(mod1, by="margin", step=10000)


# Model 2: controlling for long, lat and water depth, use mineral composition as explanatory variables.
# subset dt to remove missing data in sedimentology profile, AND remove Rhein and Ueberlingen Samples
dt_sub <- dt %>% subset_metabarlist(., "samples", !.$samples$sediment %in% 
                                      c("BS19.26","BS19.27","BS19.28","BS19.30","BS19.45","BS19.46","BS19.55")) %>%
  subset_metabarlist(., "samples",  .$samples$long > 9.2 & .$samples$long < 9.57)


sedi_sub <- sedi %>% filter(long %in% dt_sub$samples$long)

sedi_st_sub <- scale(sedi_sub %>% select(-site), center = T, scale = T) %>% as.data.frame



# Pre-select models
mod20 <- rda(dt_sub$reads ~ 1 + Condition(long+ lat + water_depth), sedi_st_sub %>% 
               select(long, lat, water_depth, Mu8.9,Qz.Mu26.6, Qz20.9, Dol31.0, Cc29.5, Chl6.2, Fsp.Mu27.8, TOC, TS, C.N.ratio),# %>% na.omit(),
             na.action=na.exclude)

mod22 <- rda(dt_sub$reads ~ TOC + TS + C.N.ratio+Chl6.2+Cc29.5+Dol31.0 + Condition(long+ lat + water_depth), sedi_st_sub, na.action = na.exclude)
vif.cca(mod22)

# select significant variables
ordistep(mod20, scope = formula(mod22), pstep=3000, R2scop=TRUE)


# model 2 with chosen variables
mod2 <- rda(dt_sub$reads ~  TS + TOC + Cc29.5 + Chl6.2 + C.N.ratio + Dol31.0 + Condition(long + lat + water_depth), sedi_st_sub, na.action=na.exclude)

# assess collinearity of variables, update mod2 by removing colinear variables if necessary
vif.cca(mod2)

# Final model 2 for all four datasets
#PLANT: Chl6.2 + Cc29.5 + Condition(long + lat + water_depth) 
#EUK: TS + TOC + Cc29.5 + Chl6.2 + C.N.ratio + Dol31.0 + Condition(long + lat + water_depth) 
#CYA: Cc29.5 + Dol31.0 + Condition(long + lat + water_depth)
#COP: Chl6.2



anova.cca(mod2, by="margin", step=10000)
RsquareAdj(mod2)$adj.r.squared

summary(mod2)

############################ Plot Supp.Fig.4 ##############################

plot_rda <- function(mod, sc=2, title) {
  ## extract % explained by the first 2 axes
  perc <- round(100*(summary(mod1)$cont$importance[2, 1:2]), 2)
  
  ## extract scores - these are coordinates in the RDA space
  sc_si <- scores(mod, display="sites", choices=c(1,2), scaling=sc)
  sc_sp <- scores(mod, display="species", choices=c(1,2), scaling=sc)
  sc_bp <- scores(mod, display="bp", choices=c(1, 2), scaling=sc)
  
  plot.new()
  par(mar=c(4,4,3,6))
  # Set up a blank plot with scaling, axes, and labels
  plot(mod,
       scaling = sc, # set scaling type 
       type = "none", # this excludes the plotting of any points from the results
       frame = FALSE, xlim=c(-1,1),
       # label the plot (title, and axes)
       main = title,
       xlab = paste0("RDA1 (", perc[1], "%)"), 
       ylab = paste0("RDA2 (", perc[2], "%)")#,
       #xaxs="i", yaxs="i"
       )
  # add points for site scores
  points(sc_si, 
         pch = 21, # set shape 
         col = "darkgoldenrod3", # outline colour "steelblue4" for 4a, "darkgoldenrod3" for 4b
         bg = "#F2bd33", # fill colour "#F2bd33" for 4b, "steelblue" for 4a
         cex = 1) # size
  # add arrows for effects of the explanatory variables
  Arrows(0,0, # start them from (0,0)
         sc_bp[,1]/1.4, sc_bp[,2]/1.4, # end them at the score value
         col = "brown4",
         lwd = 2, arr.type="triangle", arr.width=0.1)
  #add text labels for arrows
  text(x = sc_bp[,1]*1.2, # adjust text coordinate to avoid overlap with arrow tip , 0.95 for 4a, 1.2 for 4b
       y = sc_bp[,2]*1.2,
       labels = rownames(sc_bp),
       col = "brown4",
       cex = 0.8,
       font = 2)
  
  # save plot as object and return
  p <- recordPlot()
  return(p)
}

for (dt_name in c("plants_r2_embl_98","euk_r2_silva_97","cyano_r2_silva_97","cop_r5_EmblSilvaCustom_85")) {
  
  dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_norm_lgchord", sep=""))
  # match sedi data to 'samples' table, allowing 100 m max distance mismatch
  sedi <- difference_left_join(dt$samples, sedi_dt, by=c("UTM.32T.E","UTM.32T.N"),  max_dist = 100, distance_col="dist") %>%
    mutate(UTM.32T.E=UTM.32T.E.x, UTM.32T.N=UTM.32T.N.x)
  
  # remove non-numeric and useless columns
  sedi[,c("sediment","KERN","VCoarseSand","CoarseSand","GoodnessofFit","Kies",
          "UTM.32T.N.x","UTM.32T.N.y","UTM.32T.E.x","UTM.32T.E.y","dist")] <- NULL
  
  # scale numeric variables
  sedi_st <- scale(sedi %>% select(-c(site)), center = T, scale = T) %>% as.data.frame
  dt_sub <- dt %>% subset_metabarlist(., "samples", !.$samples$sediment %in% 
                                        c("BS19.26","BS19.27","BS19.28","BS19.30","BS19.45","BS19.46","BS19.55")) %>%
    subset_metabarlist(., "samples",  .$samples$long > 9.2 & .$samples$long < 9.57)
  
  
  sedi_sub <- sedi %>% filter(long %in% dt_sub$samples$long)
  
  sedi_st_sub <- scale(sedi_sub %>% select(-site), center = T, scale = T) %>% as.data.frame
  
  
  mod1 <- rda(dt$reads ~ long + lat + water_depth, sedi_st, na.action=na.exclude)
  
  if (startsWith(dt_name,"p")) {
    a1 <- plot_rda(mod1,2,"Plant")
    mod2 <- rda(dt_sub$reads ~ Chl6.2 + Cc29.5 + Condition(long + lat + water_depth), sedi_st_sub, na.action=na.exclude)
    b1 <- plot_rda(mod2,2,"Plant")
    
  } else if (startsWith(dt_name,"e")) {
    a2 <- plot_rda(mod1,2,"Eukaryote")
    mod2 <- rda(dt_sub$reads ~ TS + TOC + Cc29.5 + Chl6.2 + C.N.ratio + Dol31.0 + Condition(long + lat + water_depth), sedi_st_sub, na.action=na.exclude)
    b2 <- plot_rda(mod2,2,"Eukaryote")
    
  } else if (startsWith(dt_name, "cy")) {
    a3 <- plot_rda(mod1,2,"Cyanobacteria")
    mod2 <- rda(dt_sub$reads ~ Cc29.5 + Dol31.0 + Condition(long + lat + water_depth), sedi_st_sub, na.action=na.exclude)
    b3 <- plot_rda(mod2,2,"Cyanobacteria")
    
  } else {
    a4 <- plot_rda(mod1,2,"Copepod")
    mod2 <- rda(dt_sub$reads ~ Chl6.2 + Condition(long + lat + water_depth), sedi_st_sub, na.action=na.exclude)
    b4 <- plot_rda(mod2,2,"Copepod")
    
  }
  
}

pdf("Supp.Fig.4a_rda_mod1.pdf", width=16, height=3.2)
plot_grid(a1, a2, a3, a4, labels = NULL, nrow=1)
dev.off()

pdf("Supp.Fig.4b_rda_mod2.pdf", width=16, height=3.2)
plot_grid(b1, b2, b3, b4, labels = NULL, nrow=1)
dev.off()
############################################################################

