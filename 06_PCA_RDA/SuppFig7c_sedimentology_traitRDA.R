
# RDA using sedimentological data as explanatory variables

#library(ade4)
library(vegan)
library(factoextra) # to vis PCA
library(metabaR)
library(cowplot)
library(dplyr)
library(fuzzyjoin) # to join by coordinates of sedimentological data using inexact match

library(ggthemes)
library(shape) # to use the Arrows function


# import data
dt <- readRDS("./00_R-repository/euk_r2_silva_97_agg_norm_lgchord")

####################################################################
# Use loaded dt, run "./10_add_trait/euk_add_trait_BTL.R" to add life mode to the EUK motus table
# To continue you need objects: data
####################################################################



dt<-data


# import sedimentological data table
sedi_dt <- read.table("./06_PCA_RDA/surfsedi_sedimentological_data.txt", header=T, sep='\t')

# match sedi data to 'samples' table
sedi <- difference_left_join(dt$samples, sedi_dt, by=c("UTM.32T.E","UTM.32T.N"),  max_dist = 100, distance_col="dist") %>%
  mutate(UTM.32T.E=UTM.32T.E.x, UTM.32T.N=UTM.32T.N.x)


# remove non-numeric and useless columns
sedi[,c("sediment","KERN","VCoarseSand","CoarseSand","GoodnessofFit","Kies",
        "UTM.32T.N.x","UTM.32T.N.y","UTM.32T.E.x","UTM.32T.E.y","dist")] <- NULL

# scale numeric variables
sedi_st <- scale(sedi %>% select(-c(site)), center = T, scale = T) %>% as.data.frame




# RDA Model 2: controlling for long, lat and water depth, use mineral composition as explanatory variables.
# subset dt to remove missing data in sedimentology profile, AND remove Rhein and Ueberlingen Samples
dt_sub <- dt %>% subset_metabarlist(., "samples", !.$samples$sediment %in% 
                                      c("BS19.26","BS19.27","BS19.28","BS19.30","BS19.45","BS19.46","BS19.55")) %>%
  subset_metabarlist(., "samples",  .$samples$long > 9.2& .$samples$long < 9.57) #.$samples$long > 9.2& .$samples$long < 9.57

# subset data to selected traits
dt_sub <- subset_metabarlist(dt_sub, "motus", dt_sub$motus$Lebensform %in% c("P","B","PB","Pa","E","T"))#


sedi_sub <- sedi %>% filter(long %in% dt_sub$samples$long)

sedi_st_sub <- scale(sedi_sub %>% select(-site), center = T, scale = T) %>% as.data.frame



# model 2 with chosen variables
mod2 <- rda(dt_sub$reads ~  TS + TOC  + Chl6.2 + C.N.ratio + Dol31.0 + Condition(long + lat + water_depth), sedi_st_sub, na.action=na.exclude)

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


trait_label <- factor(dt_sub$motus$Lebensform) #Levels: B E P Pa PB T

cols = c('purple','cadetblue1','darkolivegreen3','darkred','azure3','orange')

mod=mod2;sc=2

# make plot
## extract % explained by the first 2 axes
perc <- round(100*(summary(mod)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(mod, display="sites", choices=c(1,2), scaling=sc)
sc_sp <- scores(mod, display="species", choices=c(1,2), scaling=sc)
sc_bp <- scores(mod, display="bp", choices=c(1, 2), scaling=sc)

plot.new()
par(mar=c(3,3,2,5))
# Set up a blank plot with scaling, axes, and labels
plot(mod,
     scaling = sc, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE, xlim=c(-0.1,0.12),ylim=c(-0.1,0.1),
     # label the plot (title, and axes)
     main = "RDA scaling 2\nASVs coloured by trait",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)")#,
     #xaxs="i", yaxs="i"
)
# add arrows for effects of the explanatory variables
Arrows(0,0, # start them from (0,0)
       sc_bp[,1]/8, sc_bp[,2]/8, # end them at the score value
       col = "red",
       lwd = 1.1, arr.type="triangle", arr.width=0.01, arr.length = 0.07)
# add points for site scores
points(sc_sp, 
       pch = 21, # set shape 
       col = NULL, # outline colour "steelblue4" for 4a, "darkgoldenrod3" for 4b
       bg = alpha(cols[trait_label],0.7), # fill colour "#F2bd33" for 4b, "steelblue" for 4a
       cex = 0.6) # size
#add text labels for arrows
text(x = sc_bp[,1]/6, # adjust text coordinate to avoid overlap with arrow tip , 0.95 for 4a, 1.2 for 4b
     y = sc_bp[,2]/6,
     labels = rownames(sc_bp),
     col = "red2",
     cex = 0.8,
     font = 2)
# add legend
legend(0.09,0.08, legend=c("benthic","epiphytic/epizoic","planktonic","parasitic","planktonic/benthic","terrestrial"), #
       cex=0.7,
       pch=16,
       col = cols,
       x.intersp=0.6,
       y.intersp=0.6)

# save plot as object and return
p <- recordPlot()

# output image
png('Supp.Fig.7c_sedimentology_traitRDA.png', width = 6, height = 4, units = 'in', res = 600)
plot_grid(p)
dev.off()

pdf('Supp.Fig.7c_sedimentology_traitRDA.pdf', width = 6, height = 4)
plot_grid(p)
dev.off()















