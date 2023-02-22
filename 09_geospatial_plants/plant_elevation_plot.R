library(tidyr)
library(metabaR)
library(dplyr)
library(plyr)
library(ggplot2)

source("./00_R-repository/aggregate_to_site.R")
source("./10_add_trait/plant_merge_traits_to_metabarlist.R")

dt <- readRDS("./00_R-repository/plants_r2_embl_98_agg_2rep_motuAgg")


# load info flora distribution dataset
info0 <- read.csv("./09_geospatial_plants/infoflora_export_5x5km.csv",header=T) %>% mutate(full_name=taxon) %>%
  dplyr::select(taxon_id, full_name, everything()) %>%
  separate(taxon, c("genus","species","others")," ") %>% dplyr::select(-others) 


# restrict info flora dataset to catchment area of Lake Constance, according to Abb. 1.2-1 Das Einzugsgebiet des Bodensees
info <- info0 %>% mutate(grid_x=gridcell_id%/%1000, grid_y = gridcell_id%%1000) %>%
  filter((grid_x==720 & grid_y==285) |
           (grid_x>= 715 & grid_y== 280) |
           (grid_x>=705 & grid_y==275) |
           (grid_x>=730 & grid_y==270) |
           (grid_x>=735 & grid_y==265) |
           (grid_x>=730 & grid_y==260) |
           (grid_x>=745 & grid_y==255) |
           (grid_x>=750 & grid_y==250) |
           (grid_x>=750 & grid_y==245) |
           (grid_x>=750 & grid_y==240) |
           (grid_x>=750 & grid_y==235) |
           (grid_x>=745 & grid_y==230) |
           (grid_x>=745 & grid_y==225) |
           (grid_x>=745 & grid_y==220) |
           (grid_x>=745 & grid_y==215) |
           (grid_x>=745 & grid_y==210) |
           (grid_x>=745 & grid_x<=785 & grid_y==205) |
           (grid_x>=750 & grid_x<=785 & grid_y==200) |
           (grid_x>=735 & grid_x<=795 & grid_y==195) |
           (grid_x>=735 & grid_x<= 800 & grid_y==190) |
           (grid_x>=720 & grid_x<= 795 & grid_x != 730 & grid_y==185) |
           (grid_x>=705 & grid_x<=790 & grid_y==180) |
           (grid_x>=705 & grid_x<=790 & grid_y==175) |
           (grid_x>=700 & grid_x<=790 & grid_y==170) |
           (grid_x>=690 & grid_x<=790 & grid_y==165) |
           (grid_x>=695 & grid_x<=780 & grid_y==160) |
           (grid_x>=695 & grid_x<= 770 & grid_x != 710 & grid_x != 715 & grid_y==155) |
           (grid_x>=720 & grid_x<=770 & grid_y==150) |
           (grid_x>=730 & grid_x<=770 & grid_y==145) |
           (grid_x>=750 & grid_x<=765 & grid_y==140) |
           (grid_x>=755 & grid_x<=760 & grid_y==135))


# check how many species are the only one in the genus
# then can just use genus name to match
genera1 <- info %>% dplyr::select(taxon_id, genus, species, full_name) %>% unique() %>% group_by(genus) %>% dplyr::summarise(n=n()) %>% filter(n==1) %>% dplyr::select(genus) %>% as.list %>% unlist 

print(paste(length(genera1), "out of", nrow(info %>% dplyr::select(taxon_id, genus, species, full_name) %>% unique()), "sp. are single sp. in the genus.", sep=" "))

print(paste(length(genera1), "out of", nrow(info %>% 
                                              dplyr::select(taxon_id, genus, species, full_name) %>% 
                                              unique() %>% group_by(genus) %>% dplyr::summarise(n=n())), 
            "genera have single sp.", sep=" "))


# strict match of species name + manually selected genera*
info_clean <- info %>% dplyr::select(taxon_id, genus, species, full_name) %>% unique()

genera_selected <- c("Zannichellia","Glechoma","Humulus", "Frangula", "Viscum", "Veratrum", "Parnassia", "Acorus","Hippuris","Asarum","Actaea","Parietaria")

taxa_matched <- dt$motus %>% filter(species_name %in% paste(info_clean$genus, info_clean$species, sep=" ") | 
                                      genus_name %in% genera_selected) %>%
  filter(is.na(species_name) | species_name !="Rubus idaeus") 

# subset data to include only species that has a match in the info flora database
sub <- subset_metabarlist(dt, "motus", dt$motus$TAXID %in% taxa_matched$TAXID )

# subset infoflora database
subif <- info %>% filter((paste(info$genus, info$species, sep=" ") %in% taxa_matched$species_name) |
                           (info$genus %in% taxa_matched[is.na(taxa_matched$species_name),]$genus_name))

# change TAXID of Viscum to that of Viscum album
sub$motus[sub$motus$TAXID==3971,]$TAXID <- 3972




# load pre-saved elevation data
## code to generate the file "grid_elevations_info_flora" is: generate_grid_elevations.R
grid_elevations <- readRDS("./00_R-repository/grid_elevations_info_flora")

# add elevation to infoflora dataset
subif_elev <- left_join(subif, grid_elevations, by="gridcell_id")


# aggregate pcrs by site
sub_toSite <- aggregate_to_site(sub, stats="FoD")



# number of presence among sites (or sediment extracts); add trait
nb_presence <- data.frame(nb_presence_sedi = colSums((sub_toSite$reads>0)*1)) %>% 
  mutate(TAXID=as.integer(rownames(.))) %>% left_join(., sub_toSite$motus, by="TAXID")



# combine elevation data and number of presence data in one table
tab <- left_join(subif_elev %>% mutate(tax_name=ifelse(genus %in% genera_selected, genus, paste(genus,species, sep=" "))),
                 nb_presence %>% mutate(tax_name=ifelse(is.na(species_name), genus_name, species_name)),
                 by="tax_name") %>%
  mutate(elevation_50 = round_any(median_elevation, 50)) # round up elevation to nearest 50m


# load plant trait table
trait <- read.table("./10_add_trait/plant_traits_final.txt", header=T, sep='\t', na.strings=c("","#N/A","NA"))
tab_new <- aggregate(tab$obs_nb, 
                     by=list(tax_name=tab$tax_name,
                             elevation_50=tab$elevation_50,
                             nb_presence_sedi=tab$nb_presence_sedi), 
                     FUN=sum) %>%
  dplyr::rename(., nb_obs_elev=x) %>%
  left_join(., trait %>% 
              mutate(tax_name=species_name) %>% 
              dplyr::select(tax_name, aquatic), by="tax_name")


# manually add trait for genus-level taxa
tab_new[tab_new$tax_name %in% c("Zannichellia","Acorus","Hippuris"),]$aquatic <- "aquatic"
tab_new <- tab_new %>% mutate(aquatic=ifelse(is.na(aquatic), "others", "aquatic"))


# plot
pdf("Fig.3b_plants_elevation_update.pdf", width = 7.2, height = 4)
tab_new %>% ggplot(aes(x=nb_presence_sedi, y=elevation_50, size=nb_obs_elev/4, 
                       color=aquatic)) + 
  geom_point() + 
  scale_color_manual(values = c(alpha("deepskyblue2", 0.5), alpha("#c54b8c", 0.2))) + 
  scale_size_identity(trans="sqrt", guide="legend", breaks=c(5,25,100,250), labels=c(5,25,100,250)*4) +
  labs(size="Number of\nobservations\nin catchment area") +
  scale_y_continuous(limits = c(400, 2900),breaks = seq(400, 3000, by = 500)) +
  scale_x_continuous(limits = c(0, 27),breaks = seq(0, 100, by = 5)) +
  xlab("Number of sites in lake where DNA was detected") +
  ylab("Elevation at natural occurrence / m") +
  theme_classic() +
  theme(legend.text=element_text(size=11),
        legend.title =element_text(size=13),
        text = element_text(size=rel(4)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)))#,
#panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()


