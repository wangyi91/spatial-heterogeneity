library(metabaR)
library(dplyr)
library(reshape2)
library(data.table)

plant_merge_traits_to_metabarlist <- function(dt) {
  # load trait table
  trait_table <- read.table("./10_add_trait/plant_traits_final.txt", 
             header=T, sep='\t', na.strings=c("","#N/A","NA")) %>% `row.names<-`(.$TAXID)
  
  
  # add traits to motus table
  new_motus <- dt$motus %>% inner_join(trait_table, 
                                      by = c("TAXID","SCIENTIFIC_NAME","family_name","genus_name","species_name"))
  out <- subset_metabarlist(dt, "motus", dt$motus$TAXID %in% new_motus$TAXID)
  out$motus <- new_motus %>% arrange(factor(TAXID, levels = out$motus$TAXID)) %>% 
    as.data.frame %>% `row.names<-`(.$TAXID)
  
  # organise cultivated, aquatic, alpine into one column 
  out$motus <- out$motus %>% mutate(traits = ifelse(alpine %in% c("TRUE","T"), "alpine",
                                                      ifelse(aquatic %in% c("aquatic","T"), "aquatic",
                                                             ifelse(cultivated %in% c("TRUE","T"), "cultivated",
                                                                    "others")))) %>%
    mutate(growth_form = ifelse(PlantGrowthForm %in% c("shrub","shrub/tree","tree"), "shrubs & trees",
                                ifelse(PlantGrowthForm=="graminoid", "grass & sedges",
                                       ifelse(PlantGrowthForm=="herb","non-grass herbs", PlantGrowthForm))))
  # replace NA with "uncertain" for growthforms
  out$motus <- out$motus %>% mutate(PlantGrowthForm = ifelse(is.na(PlantGrowthForm), "uncertain", PlantGrowthForm))
  
  return(out)
  }