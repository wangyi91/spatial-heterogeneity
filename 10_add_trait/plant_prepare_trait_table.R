library(metabaR)
library(dplyr)
library(reshape2)
library(data.table)

# import data, already motu-aggregated and log-chord transformed
dt_name <- "plants_r2_embl_98"
dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_motuAgg",sep=""))

# load trait table
trait <- read.table("./10_add_trait/TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease_clean.txt", 
                    header=T, sep='\t', na.strings=c("","NA"))


# Assign growth form to reads with species name assigned
tmp <- left_join(dt$motus %>% 
                   filter(!is.na(species_name)), #%>% mutate(sequence_id=rownames(.)), 
                 trait %>% 
                   mutate(species_name=AccSpeciesName) %>% 
                   select(species_name, PlantGrowthForm, Parasitic, Aquatic, Crop, Woodiness), by="species_name")

# Assign growth form to reads with genus name assigned
## since sp. in one genus can have different growth forms, here I take unique words for the list of growth forms in this genus 
tmp2 <- left_join(dt$motus %>%
                    filter(is.na(species_name)), #%>%mutate(sequence_id=rownames(.))
                  trait %>% 
                    mutate(genus_name=Genus) %>% 
                    select(genus_name, PlantGrowthForm, Parasitic, Aquatic, Crop, Woodiness) %>% 
                    group_by(genus_name) %>% 
                    dplyr::summarise(PlantGrowthForm = paste(sort(unique(na.omit(unlist(strsplit(PlantGrowthForm, "\\s+|[[:punct:]]"))))), collapse ="/"),
                              Parasitic = paste(sort(unique(na.omit(unlist(Parasitic)))), collapse ="/"),
                              Aquatic = paste(sort(unique(na.omit(unlist(Aquatic)))), collapse ="/"),
                              Crop = paste(sort(unique(na.omit((Crop)))), collapse ="/"),
                              Woodiness = paste(sort(unique(na.omit(unlist(strsplit(Woodiness,"/"))))), collapse ="/")), #unlist(strsplit(na.omit(PlantGrowthForm), "\\s+|[[:punct:]]"))
                  by="genus_name")

motus_trait <- bind_rows(tmp, tmp2)


# replace "" with NA
motus_trait <- motus_trait %>% mutate_all(na_if,"")

# update the motus table
dt$motus <- left_join(dt$motus, #%>% mutate(sequence_id=rownames(.))
                      motus_trait %>% select(TAXID, PlantGrowthForm, Parasitic, Aquatic, Crop, Woodiness), 
                      by="TAXID")

rownames(dt$motus) <- dt$motus$TAXID
identical(rownames(dt$motus), colnames(dt$reads))
check_metabarlist(dt)

# subset data to keep taxa identified to species level
sub <- subset_metabarlist(dt, "motus", !is.na(dt$motus$species_name))

# Manually mark alpine taxa
alp_ls <- c("Androsace","Anthyllis montana","Arctostaphylos uva-ursi", "Arctous alpina",     
"Bartsia alpina","Calluna vulgaris","Cerastium cerastoides","Chrysosplenium wrightii",
"Doronicum","Edraianthus","Kalmia procumbens","Myricaria germanica","Oxyria digyna",
"Pedicularis","Phyllodoce nipponica","Phyteuma","Ranunculus glacialis",
"Ranunculus reptans","Saxifraga aspera","Saxifraga berica","Saxifraga callosa","Saxifraga kotschyi","Saxifraga lingulata","Saxifraga paniculata","Saxifraga taygetea","Senecio",
"Trachydium variabile","Trifolium pallescens","Vaccinium vitis-idaea","Valeriana elongata",
"Veronica","Veronica alpina","Veronica aphylla","Viola brevistipulata") 

# load data from info flora (Flora Helvetica)
IF <- read.csv("./10_add_trait/info_flora_Checklist_simplified.csv", header=T, na.strings="") %>% filter(Rangstufe == "sp")


# mark alpine plants &
# mark what species are present in info flora checklist
sub$motus <- sub$motus %>% mutate(alpine = sub$motus$genus_name %like% paste(alp_ls, collapse  = "|") | 
                                    sub$motus$species_name %like% paste(alp_ls, collapse  = "|"),
                                  info_flora = sub$motus$species_name %in% paste(IF$Gattung, IF$Spezies, sep=" "))


# create a column marking cultivated or non-cultivated plants in infoflora
## get a species list
IF_sp <- sub$motus[sub$motus$info_flora,]$species_name %>% unique

## non-cultivated plants: A, I, ni
ncul <- IF %>% filter(paste(IF$Gattung, IF$Spezies, sep=" ") %in% IF_sp) %>% filter(Indigenat.CH %in% c("I","A","I/N","A/N","ni"))

## cultivated plants: AC, NC
neo_cul <- IF %>% filter(paste(IF$Gattung, IF$Spezies, sep=" ") %in% IF_sp) %>% filter(Indigenat.CH %in% c(NA,"NC","AC","N"))

sub$motus <- sub$motus %>%
  mutate(cultivated_IF =  ifelse(sub$motus$species_name %in% paste(neo_cul$Gattung, neo_cul$Spezies, sep=" "), "TRUE", 
                                 ifelse(sub$motus$species_name %in% paste(ncul$Gattung, ncul$Spezies, sep=" "), "FALSE", NA)))
                                  

# manually supplement cultivated species, aquatic species and missing growth forms:
cult_ls <- c("Acer pictum","Acer pycnanthum","Acer saccharum",        
"Aesculus pavia","Ananas comosus","Arctotis breviscapa",
"Artemisia armeniaca","Arundinaria fargesii","Asimina triloba",
"Avena macrostachya","Barclaya longifolia","Berberis aristata",
"Berberis japonica","Calocedrus decurrens","Cannabis sativa",
"Cedrus deodara","Chamaenerion angustifolium", "Citrus sinensis",
"Corchorus olitorius","Cornus florida","Curtisia dentata",
"Davidia involucrata","Dendropanax hainanensis","Dryopteris erythrosora",
"Ficus insipida","Fragaria viridis","Ginkgo biloba",
"Indigofera tinctoria","Juniperus chinensis","Lonicera involucrata",
"Lyonothamnus floribundus","Nerium oleander","Passiflora incarnata",
"Pemphis acidula","Picea breweriana","Picea neoveitchii",
"Platanus occidentalis","Platanus racemosa","Poa ligulata",
"Prunus oblonga","Quercus mongolica","Ricinus communis",
"Rubus occidentalis","Suriana maritima","Tecoma stans",
"Theobroma cacao","Thujopsis dolabrata","Tilia mandshurica",
"Viola brevistipulata")

aqua_ls <- c("Stuckenia pectinata", "Stuckenia filiformis","Potamogeton octandrus")

growthform_ls <- read.table("./10_add_trait/growthform_manual_list.txt", header=T, sep='\t', na.strings="NA")

sub$motus <- sub$motus %>%
  mutate(cultivated = ifelse(species_name %in% cult_ls, "T", cultivated_IF),
         aquatic = ifelse(sequence %in% aqua_ls, "aquatic", Aquatic),
         PlantGrowthForm = ifelse(species_name %in% growthform_ls$species_name, growthform_ls$growth_form, PlantGrowthForm))  





trait_out <- sub$motus %>% select(TAXID,SCIENTIFIC_NAME,family_name,genus_name,species_name,
                                        PlantGrowthForm,aquatic,alpine,info_flora,cultivated) #Parasitic,Crop,Woodiness,


# out put motus table to sort traits manually
write.table(trait_out, './10_add_trait/plant_traits_final.txt', row.names=FALSE, col.names=T, sep = "\t", quote = FALSE)
