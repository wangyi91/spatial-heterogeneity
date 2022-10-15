# RDA using longitude, latitude and sedimentological data as explanatory variables

library(vegan) # function rda
library(metabaR)
library(dplyr)

#import data
dt_clean<-readRDS("./00_R-repository/plants_r2_embl_98_clean")
dt_clean<-readRDS("./00_R-repository/euk_r2_silva_97_agg")
dt_clean<-readRDS("./00_R-repository/cyano_r2_silva_97_agg")
dt_clean<-readRDS("./00_R-repository/cop_r2_EmblSilvaCustom_85_agg")




# aggregate pcrs by core, output detection probability
dt_cr <- dt_clean
dt_cr$pcrs <- left_join(dt_clean$pcrs, dt_clean$samples %>% mutate(sample_id=rownames(.)), by=c("sediment","sample_id"))
rownames(dt_cr$pcrs) <- rownames(dt_clean$pcrs)
dt_cr$pcrs$sample_id <- dt_cr$pcrs$sediment; dt_cr$pcrs$sediment <- NULL

cores<- read.table("./03_metabaR/surfsedi_cores.txt", header=T, sep='\t')
rownames(cores) <- cores$sediment
dt_cr$samples <- cores

dt <- aggregate_pcrs(dt_cr, FUN=FUN_agg_pcrs_prob) #%>%
#subset_metabarlist(., "samples",  .$samples$long > 9.2)# & dt$samples$long<9.57) 
# exclude Rhein inflow samples (outlier) and ueberlingen samples (no sedimentological data)






# import sedimentological data table
sedi_dt <- read.table("./06_redundancy_analysis/surfsedi_sedimentological_data.txt", header=T, sep='\t')

# match sedi data to 'samples' table, allowing 100 m max distance mismatch
sedi <- difference_left_join(dt$samples, sedi_dt, by=c("UTM.32T.E","UTM.32T.N"),  max_dist = 100, distance_col="dist") %>%
  mutate(UTM.32T.E=UTM.32T.E.x, UTM.32T.N=UTM.32T.N.x)

rownames(sedi) <- sedi$sediment

# remove non-numeric and useless columns
sedi[,c("sediment","KERN","VCoarseSand","CoarseSand","GoodnessofFit","Kies",
        "UTM.32T.N.x","UTM.32T.N.y","UTM.32T.E.x","UTM.32T.E.y","dist")] <- NULL





# check collinearity between variables
sub <- sedi %>% 
  select(long, lat, water_depth, Sand,Mud,Clay,Mu8.9:Dol31.0,TS,TOC,TN) %>% 
  arrange(long) %>% 
  mutate(id=rownames(.)) #%>% select(-c(   Chl6.2, Qz.Mu26.6, Fsp.Mu27.8))

sub_long <- reshape2::melt(sub, id.vars="id") %>% na.omit() #%>% filter(variable %in% c("TS"  )) #"Mu8.9" ,"Chl6.2","Qz20.9", "Fsp.Mu27.8", "Dol31.0"  "Qz.Mu26.6",  "Cc29.5", 

ggplot(sub_long, aes(x=id,y=value,col=variable,group=variable)) + geom_line() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
















