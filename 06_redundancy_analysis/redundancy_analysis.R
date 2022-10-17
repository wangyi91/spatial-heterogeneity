# RDA using longitude, latitude and sedimentological data as explanatory variables

library(vegan) # function rda
library(metabaR)
library(dplyr)
library(fuzzyjoin)

# import data. Choose one
dt_name = "plants_r2_embl_98"
dt_name = "euk_r2_silva_97"
dt_name = "cyano_r2_silva_97"
dt_name = "cop_r5_EmblSilvaCustom_85"

dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_norm_lgchord", sep=""))



# import sedimentological data table
sedi_dt <- read.table("./06_redundancy_analysis/surfsedi_sedimentological_data.txt", header=T, sep='\t')

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



# Model 2: controllong for long, lat and water depth, use mineral composition as explanatory variables.
# subset dt to remove missing data in sedimentology profile, AND remove Rhein and Ueberlingen Samples
dt_sub <- dt %>% subset_metabarlist(., "samples", !.$samples$sediment %in% 
                                      c("BS19.26","BS19.27","BS19.28","BS19.30","BS19.45","BS19.46","BS19.55")) %>%
  subset_metabarlist(., "samples",  .$samples$long > 9.2 & .$samples$long < 9.57)


sedi_sub <- sedi %>% filter(long %in% dt_sub$samples$long)

sedi_st_sub <- scale(sedi_sub %>% select(-site), center = T, scale = T) %>% as.data.frame



# Pre-select models
mod20 <- rda(dt_sub$reads ~ 1 + Condition(long+ lat + water_depth), sedi_st_sub %>% 
               select(long, lat, water_depth, Mu8.9,Qz.Mu26.6, Qz20.9, Dol31.0, Cc29.5, Chl6.2, Fsp.Mu27.8, TOC, TS, C.N.ratio) %>% na.omit(),
             na.action=na.exclude)

mod22 <- rda(dt_sub$reads ~ TOC + TS + C.N.ratio+Chl6.2+Cc29.5+Dol31.0 + Condition(long+ lat + water_depth), sedi_st_sub %>% na.omit, na.action = na.exclude)
vif.cca(mod22)

# select significant variables
ordistep(mod20, scope = formula(mod22), pstep=3000, R2scop=TRUE)


# model 2 with chosen variables
mod2 <- rda(dt_sub$reads ~  Chl6.2  + Condition(long + lat + water_depth), sedi_st_sub, na.action=na.exclude)

# assess collinearity of variables, update mod2 by removing colinear variables if necessary
vif.cca(mod2)

# Final model 2 for all four datasets
#PLANT: Chl6.2 + Cc29.5 + Condition(long + lat + water_depth) 
#EUK: TS + TOC + Cc29.5 + Chl6.2 + C.N.ratio + Dol31.0 + Condition(long + lat + water_depth) 
#CYA: Cc29.5 + Dol31.0 + Condition(long + lat + water_depth)
#COP: Chl6.2



anova.cca(mod2, by="margin", step=5000)
RsquareAdj(mod2)$adj.r.squared

summary(mod2)





