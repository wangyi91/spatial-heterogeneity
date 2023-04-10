
# RDA using sedimentological data as explanatory variables

library(vegan)
library(factoextra) # to vis PCA
library(metabaR)
library(cowplot)
library(dplyr)
library(fuzzyjoin) # to join by coordinates of sedimentological data using inexact match

library(ggthemes)
library(shape) # to use the Arrows function

source("./00_R-repository/aggregate_to_site.R")

# import data
dt <- readRDS("./00_R-repository/euk_r2_silva_97_agg_norm_lgchord")
#dt <- readRDS("./00_R-repository/euk_r2_silva_97_agg_2rep_motuAgg_norm")

# rate of detection at sites
#dt <- aggregate_to_site(dt_clean, stats="FoD")

#subset_metabarlist(., "samples",  .$samples$long > 9.2)# & dt$samples$long<9.57) 
# exclude Rhein inflow samples (outlier) and ueberlingen samples (no sedimentological data)


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




######################## Plot Author_Response_Fig3a ####################################

# sites to plot sedimentological data. Sites (S13, 15, 18) with NA are excluded
ss <- c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11","S12","S14","S16",
        "S17","S19","S20","S22") 



psite <- sedi %>% select(-c(water_depth, long, lat, UTM.32T.E,UTM.32T.E.dist,UTM.32T.N,UTM.32T.N.dist)) %>%
  mutate(Silt=Mud-Clay) %>% 
  filter(site %in% ss) %>% distinct() 
psite_long <- reshape2::melt(psite, id.vars="site") %>% na.omit() #%>% filter(variable %in% c("Chl6.2","Fsp.Mu27.8","Dol31.0","Cc29.5","Qz.Mu26.6", "TS","TOC","TN"))

plot1 <- psite_long %>% filter(variable %in% c("Chl6.2","Fsp.Mu27.8","Dol31.0","Cc29.5","Qz.Mu26.6")) %>%
  ggplot(aes(x=factor(site, levels=rev(ss)), y=value, col=variable, 
             group=factor(variable, levels = c("Chl6.2","Fsp.Mu27.8","Dol31.0","Cc29.5","Qz.Mu26.6")))) +#"Chl6.2","Qz20.9","Dol31.0","Cc29.5","Qz.Mu26.6", 
  geom_col(colour="black", size=0.2, width=0.75, aes(fill=variable)) + 
  xlab("site") + 
  ylab("count by X-ray diffraction, k/s") + 
  labs(fill='Minerals') +
  scale_fill_manual(labels=c("Chlorite", 
                             "Quartz-muscovite",
                             "Calcite",
                             "Feldspar-muscovite",
                             "Dolomite"), values=rev(c("#90005D","#E5A9B4", "#F8FCE3", "#B5BE7E", "#798233")))+ #hcl.colors(5,"ArmyRose"))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))



plot2 <- psite_long %>% filter(variable %in% c("TS","TOC","TN")) %>%
  ggplot(aes(x=factor(site, levels=rev(ss)), y=value, 
             col=variable, 
             group=factor(variable, levels = c("TS","TN","TOC")))) + 
  geom_col(colour="black", size=0.2, width=0.75, aes(fill=variable)) + 
  ylab("%") + 
  labs(fill='Elements') +
  scale_fill_manual(labels=c("Total nitrogen", 
                             "Total organic carbon",
                             "Total sulphur"), values=rev(hcl.colors(3,"Purple-Blue"))) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_reverse() +
  scale_x_discrete(position="top")

plot3 <- psite_long %>% 
  filter(variable %in% c("Sand","Silt","Clay")) %>%
  ggplot(aes(x=factor(site, levels=rev(ss)), y=value, 
             col=variable, 
             group=factor(variable, levels = c("Sand","Silt","Clay")))) + 
  geom_col(colour="black", size=0.2, width=0.75, aes(fill=variable)) + 
  ylab("%") + 
  labs(fill='Particles') +
  scale_fill_manual(labels=c("Sand","Silt","Clay"), values=rev(hcl.colors(3,"Fall")[c(2,1,3)])) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_reverse() +
  scale_x_discrete(position="top")


# align the two plots
cowplot::plot_grid(plot3, plot2, plot1, align = "v", ncol = 1, rel_heights = c(0.25, 0.25,0.5))


# output image
png('Supp.Fig.7a_sedimentology.png', width = 7.2, height = 4.8, units = 'in', res = 600)
cowplot::plot_grid(plot3, plot2, plot1, align = "v", ncol = 1, rel_heights = c(0.25, 0.25,0.5))
dev.off()

pdf('Supp.Fig.7a_sedimentology.pdf', width = 7.2, height = 4.8)
cowplot::plot_grid(plot3, plot2, plot1, align = "v", ncol = 1, rel_heights = c(0.25, 0.25,0.5))
dev.off()

##########################################################################################


# RDA Model 2: controlling for long, lat and water depth, use mineral composition as explanatory variables.
# subset dt to remove missing data in sedimentology profile, AND remove Rhein and Ueberlingen Samples
dt_sub <- dt %>% subset_metabarlist(., "samples", !.$samples$sediment %in% 
                                      c("BS19.26","BS19.27","BS19.28","BS19.30","BS19.45","BS19.46","BS19.55")) %>%
  subset_metabarlist(., "samples",  .$samples$long > 9.2& .$samples$long < 9.57) #.$samples$long > 9.2& .$samples$long < 9.57


# subset data to selected taxa
dt_sub <- dt_sub %>% subset_metabarlist(., "motus", (!is.na(.$motus$phylum_name)&
                               .$motus$phylum_name %in% 
                                 c("Annelida","Nematoda","Arthropoda", "Platyhelminthes","Bacillariophyta","Chlorophyta")) |
                                 (!is.na(.$motus$kingdom_name)&.$motus$kingdom_name =="Fungi"))#


sedi_sub <- sedi %>% filter(long %in% dt_sub$samples$long)

sedi_st_sub <- scale(sedi_sub %>% select(-site), center = T, scale = T) %>% as.data.frame



# Pre-select models
mod20 <- rda(dt_sub$reads ~ 1 + Condition(long + lat + water_depth), 
             sedi_st_sub %>% 
               select(long, lat, water_depth, Mu8.9,Qz.Mu26.6, Qz20.9, Cc29.5, Chl6.2, Dol31.0,Fsp.Mu27.8, TOC, TS, C.N.ratio) %>% na.omit(),
             na.action=na.exclude) 

mod22 <- rda(dt_sub$reads ~ TOC + TS + C.N.ratio+Chl6.2+Cc29.5 + Dol31.0 + Condition(long+ lat + water_depth), sedi_st_sub, na.action = na.exclude)
vif.cca(mod22)

# select significant variables
ordistep(mod20, scope = formula(mod22), pstep=3000, R2scop=TRUE)


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

tax_label <- factor(ifelse((!is.na(dt_sub$motus$kingdom_name)&dt_sub$motus$kingdom_name=="Fungi"), "Fungi", dt_sub$motus$phylum_name))
# Levels: Annelida Arthropoda Bacillariophyta Chlorophyta Fungi Nematoda Platyhelminthes
cols = c('darkred','gold','cadetblue1','darkolivegreen3','purple','azure3','orange')

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
     frame = FALSE, xlim=c(-0.08,0.12),ylim=c(-0.08,0.08),
     # label the plot (title, and axes)
     main = "RDA scaling 2\nASVs coloured by taxon",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)")#,
     #xaxs="i", yaxs="i"
)
# add points for site scores
points(sc_sp, 
       pch = 21, # set shape 
       col = NULL, # outline colour "steelblue4" for 4a, "darkgoldenrod3" for 4b
       bg = alpha(cols[tax_label],0.7), # fill colour "#F2bd33" for 4b, "steelblue" for 4a
       cex = 0.8) # size
# add arrows for effects of the explanatory variables
Arrows(0,0, # start them from (0,0)
       sc_bp[,1]/8, sc_bp[,2]/8, # end them at the score value
       col = "red",
       lwd = 2, arr.type="triangle", arr.width=0.03, arr.length = 0.1)
#add text labels for arrows
text(x = sc_bp[,1]/6, # adjust text coordinate to avoid overlap with arrow tip , 0.95 for 4a, 1.2 for 4b
     y = sc_bp[,2]/6,
     labels = rownames(sc_bp),
     col = "red2",
     cex = 0.8,
     font = 2)
# add legend
legend(0.09,0.08, legend=levels(tax_label), 
       cex=0.7,
       pch=16,
       col = cols,
       x.intersp=0.6,
       y.intersp=0.6)

# save plot as object and return
p <- recordPlot()

# output image
png('Supp.Fig.7b_sedimentology_taxRDA.png', width = 6, height = 4, units = 'in', res = 600)
plot_grid(p)
dev.off()

pdf('Supp.Fig.7b_sedimentology_taxRDA.pdf', width = 6, height = 4)
plot_grid(p)
dev.off()















