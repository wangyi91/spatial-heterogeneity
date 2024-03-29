library(metabaR)
library(dplyr)
library(stringr)
library(data.table)
library(readxl)

#import data, already motu-aggregated and log-chord transformed
dt <- readRDS("./00_R-repository/euk_r2_silva_97_agg_2rep_norm")

dt$motus$sequence_id <- rownames(dt$motus)


# load trait table, only keep lifeform column
trait_full <- read_xlsx("./10_add_trait/BTL_Mai_2020.xlsx") %>% 
  `colnames<-`(sub(" ", "_", colnames(.))) %>% `colnames<-`(sub("-", "_", colnames(.))) %>%
  mutate(Taxon=ifelse(Status=="Synonym", gültiges_Taxon, Taxon))

# if a sp has more than one possible living form, exclude for now
trait <- trait_full %>% select(Taxon, Lebensform_Typ) %>% filter(!is.na(Lebensform_Typ)) %>% unique %>% 
  group_by(Taxon) %>% mutate(n=n()) %>% ungroup %>%
  filter(n<2 | (Taxon %in% c("	
Navicula cryptocephala","Ulnaria ulna") & Lebensform_Typ !="PB")) %>% select(-n)


# output taxon names to get full lineage using TaxonKit
# 1. full species name
# 2. if 1. not available, get genus name i.e. sudo_genus
taxa_search <- trait %>% mutate(taxon_to_search = ifelse(is.na(word(trait$Taxon,1,2)), word(trait$Taxon,1), word(trait$Taxon,1,2))) %>% select(-Taxon) %>% unique

write.table(taxa_search$taxon_to_search, "./taxon_names.txt",quote=F,sep="\t", row.names = F, col.names = F)

# Run TaxonKit
system("cd [YOUR_PATH]/03_codes/17_add_trait/")

# run TaxonKit
system(
  "cat taxon_names.txt | taxonkit name2taxid | taxonkit lineage -i 3 -r > lineage.txt"
)


# then manually delete rows in lineage.txt where there is no taxid
lineage <- read.csv("./lineage.txt", sep="\t", header=F) %>% filter(!is.na(V2))
write.table(lineage, "./lineage_rev.txt", quote=F, sep="\t", col.names = F, row.names = F)


# use Taxonkit to add family and phylum info
system(
  'cat lineage_rev.txt | taxonkit reformat -I 2 -r NA -f "{p}\t{f}\t{g}" > pfg.txt'
)


f <- read.table("./10_add_trait/pfg.txt", header=F,sep="\t", na.strings = c("NA",""),
                col.names = c("taxon_to_search","TAXID","path","rank","phylum_name","family_name","genus_name"))

trait_join <- left_join(trait %>% mutate(taxon_to_search = 
                                           ifelse(is.na(word(trait$Taxon,1,2)), 
                                                  word(trait$Taxon,1), 
                                                  word(trait$Taxon,1,2))) , 
                        f %>% select(taxon_to_search,phylum_name, family_name, genus_name) %>% unique,
                        by="taxon_to_search")

trait_auto1 <- trait_join %>%
  filter(!is.na(phylum_name)|!is.na(family_name)|!is.na(genus_name))


# for sp. in trait_auto that have no match from Taxonkit search, change taxon_to_search to genus i.e. sudo_genus and search again
taxa_search2 <- trait_join %>% filter(is.na(phylum_name), is.na(family_name), is.na(genus_name)) %>% mutate(taxon_to_search=word(taxon_to_search,1)) %>% select(taxon_to_search) %>% unique

write.table(taxa_search2$taxon_to_search, "./taxon_names_2.txt",quote=F,sep="\t", row.names = F, col.names = F)


# run TaxonKit
system(
  "cat taxon_names_2.txt | taxonkit name2taxid | taxonkit lineage -i 3 -r > lineage2.txt")


# then manually delete rows in lineage.txt where there is no taxid
lineage2 <- read.csv("./lineage2.txt", sep="\t", header=F) %>% filter(!is.na(V2))
write.table(lineage2, "./lineage2_rev.txt", quote=F, sep="\t", col.names = F, row.names = F)


# use Taxonkit to add family and phylum info
system(
  'cat lineage2_rev.txt | taxonkit reformat -I 2 -r NA -f "{p}\t{f}\t{g}" > pfg2.txt')


f2 <- read.table("./10_add_trait/pfg2.txt", header=F,sep="\t", na.strings = c("NA",""),
                 col.names = c("taxon_to_search","TAXID","path","rank","phylum_name","family_name","genus_name")) %>%
  filter(rank=="genus") 

trait_auto2 <- left_join(trait %>% filter(!Taxon %in% trait_auto1$Taxon) %>% 
                           mutate(taxon_to_search = word(.$Taxon,1)), 
                         f2 %>% select(taxon_to_search,phylum_name, family_name, genus_name) %>% unique,
                         by="taxon_to_search")


# combine auto1 and 2
trait_auto <- bind_rows(trait_auto1,trait_auto2)




# load manually-added life form trait table
trait_man <- read_xlsx("./10_add_trait/euk_trait_manual.xlsx", na = "NA", sheet=1) %>% 
  filter(!is.na(Lebensform_manual) & !is.na(family_name)) %>%
  select(phylum_name, family_name, genus_name, SCIENTIFIC_NAME, Lebensform_manual)


# combine trait_auto and trait_man
trait_comb <- bind_rows(trait_auto %>% rename(Lebensform = Lebensform_Typ),
                        trait_man %>% rename(Taxon=SCIENTIFIC_NAME, Lebensform = Lebensform_manual)) %>% 
  distinct(Lebensform, taxon_to_search, genus_name, family_name, phylum_name, .keep_all = TRUE)


# Assign growth form to reads with species name assigned
## remove poorly assigned species such as 'environmental sample', 'sp.' and 'cf.', and those with non-latin words and those without genus name.
motus_clean <- dt$motus %>% 
  mutate(poor_sp_assign=as.integer(species_name %like% 'sample' | species_name %like% 'sp.'| 
                                     species_name %like% 'cf.' |
                                     species_name %like% 'uncultured' | species_name %like% 'Spumella-like' |
                                     species_name %like% 'green ' | is.na(genus_name))) %>%
  mutate(species_name=ifelse(poor_sp_assign==1, NA, species_name)) %>%
  mutate(sequence_id=rownames(.))



tmp1 <- left_join(motus_clean %>% filter(!is.na(species_name) & !is.na(genus_name)), 
                  trait_comb %>% mutate(species_name=word(Taxon,1,2)) %>% 
                    select(-family_name, -phylum_name,-genus_name, -Taxon, -taxon_to_search) %>% unique, 
                  by="species_name") %>% 
  group_by(sequence) %>% mutate(n=n()) %>% filter(!(n>1&Lebensform=="B")) %>% # this is done after visual check
  ungroup %>% select(-n)


# Assign growth form to reads with genus name assigned
## since sp. in one genus can have different growth forms, here I take unique words for the list of living forms in this genus, and keep those with only one possible value
tmp2 <- left_join(motus_clean %>% filter(is.na(species_name) & !is.na(genus_name)), 
                  trait_comb %>% group_by(genus_name) %>%
                    summarise(Lebensform=unique(Lebensform)) %>% 
                    filter(!is.na(Lebensform)) %>% mutate(n=n()) %>% ungroup %>% filter(n==1)%>% select(-n), 
                  by="genus_name")

# family 
tmp3 <- left_join(motus_clean %>% filter(is.na(genus_name) & !is.na(family_name)), 
                  trait_comb %>% group_by(family_name) %>%
                    summarise(Lebensform=unique(Lebensform)) %>% 
                    filter(!is.na(Lebensform)) %>% mutate(n=n()) %>% ungroup %>% filter(n==1)%>% select(-n), 
                  by="family_name")

motus_trait <- bind_rows(tmp1, tmp2, tmp3) %>% filter(!is.na(Lebensform)) %>% as.data.frame
rownames(motus_trait) <- motus_trait$sequence_id





## subset dt to include only ASV with trait assigned
data <- subset_metabarlist(dt, "motus", rownames(dt$motus) %in% motus_trait$sequence_id)

data$motus <- left_join(data$motus, 
                        motus_trait %>% select(sequence_id, Lebensform), 
                        by="sequence_id")

rownames(data$motus) <- data$motus$sequence_id

check_metabarlist(data)



