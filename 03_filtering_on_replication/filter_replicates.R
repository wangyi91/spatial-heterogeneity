# This script filters out ASVs that appear in too few replicates. 
# Then it aggregate ASVs of the same taxonomic identity.

library(metabaR)
library(dplyr)

#import data
dt_name <- "plants_r2_embl_98"
dt_name <- "euk_r2_silva_97"
dt_name <- "cyano_r2_silva_97"
dt_name <- "cop_r5_EmblSilvaCustom_85"

dt <- readRDS(paste("./00_R-repository/",dt_name,"_clean", sep=""))

dt_agg <- aggregate_pcrs(dt, FUN=FUN_agg_pcrs_mean) 

#dt_agg_motuAgg <- aggregate_motus(dt_agg, dt_agg$motus$TAXID)

saveRDS(dt_agg, file = paste("./00_R-repository/", dt_name, "_agg", sep=""))

#saveRDS(dt_agg_motuAgg, file = paste("./00_R-repository/", dt_name, "_agg_motuAgg", sep=""))



# criteria to keep ASVs: appear at least in r/3 replicates
r=2

# get reads who are present in more than 1/3 PCR replicates
dt_pa <- dt
dt_pa$reads <- (dt$reads>0)*1 
dt_pa_agg <- aggregate_pcrs(dt_pa, FUN=FUN_agg_pcrs_prob)

# mark ASVs with >= r replicates as 1, else 0
index_nRep <- (dt_pa_agg$reads>=r/3)*1

dt_agg_nrep <- dt_agg
dt_agg_nrep$reads <- dt_agg$reads * index_nRep

# clean up metabarlist
motus_keep <- colSums(dt_agg_nrep$reads) %>% as.data.frame %>% filter(.>0) %>% rownames
dt_agg_nrep <- subset_metabarlist(dt_agg_nrep, "motus", rownames(dt_agg_nrep$motus) %in% motus_keep)
check_metabarlist(dt_agg_nrep)

saveRDS(dt_agg_nrep, file = paste("./00_R-repository/", dt_name, "_agg_", r, "rep", sep=""))

# aggregate motus by TAXID and save file
dt_motuAgg_nrep <- aggregate_motus(dt_agg_nrep, groups = dt_agg_nrep$motus$TAXID)
saveRDS(dt_motuAgg_nrep, file = paste("./00_R-repository/", dt_name, "_agg_",r,"rep_motuAgg", sep=""))


