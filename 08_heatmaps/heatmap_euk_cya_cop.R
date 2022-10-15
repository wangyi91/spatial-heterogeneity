# Create heatmaps for taxonomic groups in the EUK, CYA and COP datasets.

library(metabaR)
library(dplyr)
library(stringr) # for string manipulation
library(gplots) # to use heatmap.2
library(RColorBrewer) # to visualise color palette
suppressPackageStartupMessages(library(dendextend)) # to manage clustering dendrogram
library(colorspace)# to work with dendextend

library(ape)
library(tidytree)
library(tibble)
library(phylogram)

library(seqRFLP)

source("./00_R-repository/aggregate_motus_geneious.R")
source("./00_R-repository/make_new_tree_geneious.R")

# import data. Choose one
dt_name = "euk_r2_silva_97"
dt_name = "cyano_r2_silva_97"
dt_name = "cop_r5_EmblSilvaCustom_85"

dt <- readRDS(paste("./00_R-repository/", dt_name, "_clean", sep=""))



# add sequence_id into motus table, for later use
dt$motus <- dt$motus %>% mutate(sequence_id=rownames(.)) %>% select(sequence_id, everything())



# choose taxonomic group of interest
tax_name="Platyhelminthes";rank_name = "phylum_name"
tax_name="Annelida";rank_name = "phylum_name"
tax_name="Nematoda";rank_name = "phylum_name"
tax_name="Bacillariophyta";rank_name = "phylum_name"
tax_name="Fungi";rank_name = "kingdom_name"
tax_name="Chlorophyta";rank_name = "phylum_name"

tax_name="Cyanobacteria";rank_name = "phylum_name"
tax_name="Copepoda";rank_name = "subclass_name"

# for Fig.S3
tax_name="Ciliophora";rank_name = "phylum_name"
tax_name="Cercozoa";rank_name = "phylum_name"
tax_name="Arthropoda";rank_name = "phylum_name"




# subset data
dt_sub <- subset_metabarlist(dt, "motus", 
                             dt$motus[,rank_name] == tax_name & !is.na(dt$motus[,rank_name]))




# filter data based on identification of chosen rank
# here I choose class as the rank to filter the EUK dataset, family for CYA dataset 
filter_rank <- ifelse(tax_name=="Cyanobacteria","family_name", "class_name")

dt_sub <- subset_metabarlist(dt_sub, "motus", 
                             !is.na(dt_sub$motus[,filter_rank]))

# use geneious prime to generate MOTUs
dt_geneious <- aggregate_motus_geneious(dt_sub, tax_name)


# aggregate pcrs by site, output frequency of detection 
dt_st <- dt_geneious
dt_st$pcrs <- left_join(dt_geneious$pcrs, dt_geneious$samples %>% mutate(sample_id=rownames(.)), by="sample_id")
rownames(dt_st$pcrs) <- rownames(dt_geneious$pcrs)
dt_st$pcrs$sample_id <- dt_st$pcrs$site

# load the new 'samples' table: sites
sites<- read.table("./03_metabaR/surfsedi_sites.txt", header=T, sep='\t')
rownames(sites) <- sites$site
dt_st$samples <- sites

# frequency of detection at sites
data <- aggregate_pcrs(dt_st, FUN=FUN_agg_pcrs_prob)




ab <- data$reads # frequency of detection

# clustering method for hclust
meth <- "complete"


# ordinary heatmap

ht <- ab %>% heatmap.2(hclustfun=function(x) hclust(x, method = meth), 
            labCol = NA, labRow = NA, trace="none", 
            col=hcl.colors(12, "Oslo", rev = F)) # add Rowv=NA to stop row from clustering


# save the order in dendrogram for later use
ddorder <- ht$colInd


# Plot the heatmap where ASVs are organised in phylogenetic tree.
# First make the new tree where tips were collapsed
make_new_tree_geneious(data, tax_name)


# rooted new tree
nx_new_rt <- read.nexus(file = paste(tax_name, "_new.nex",sep = "")) # :- replaced to _






# Make taxonomy label of chosen rank to plot in heatmap
# choose rank to plot in heatmap
if (tax_name=="Cyanobacteria"){
  rank_label <- "family_name"
} else if (tax_name=="Copepoda"){
  rank_label <- "genus_name"
} else {
  rank_label <- "order_name"
}


# match nx_new_rt$tip.label to colnames of "ab", the former of which heatmap can use in plotting tree in heatmap
new_tip <- left_join(nx_new_rt$tip.label %>% as.data.frame %>% `colnames<-`("label"), 
                     data$motus %>% mutate(label=gsub(":","_",substring(sequence_id,29))) %>% 
                       mutate(id=rownames(.)) %>% select(label,id, all_of(rank_label), BEST_IDENTITY),
                     by="label")

new_tip[!is.na(new_tip[,"id"]) & is.na(new_tip[,rank_label]) | (!is.na(new_tip[,"BEST_IDENTITY"]) & new_tip[,"BEST_IDENTITY"] < 0.85),rank_label]="(unidentified)"


# further sort reads using tax information, when tree and identification are not congruent.
if (dt_name %in% c("euk_r2_silva_97","cyano_r2_silva_97")) {
  motus_sorted <- data$motus %>% dplyr::slice(all_of(ddorder)) %>% #first order rows as in dendrogram so that it is kept later
    select(sequence_id, class_name, order_name, family_name, genus_name) %>%
    arrange(., is.na(get(rank_label)), get(filter_rank), get(rank_label)) %>%
    mutate(label=gsub(":","_",substring(sequence_id,29)))
  
  new_tip <- right_join(motus_sorted %>% select(-all_of(rank_label)), new_tip , by= "label")
} else {
  # if working on COP:
  new_tip <- right_join(data$motus %>% slice(all_of(ddorder)) %>%
                        mutate(label=gsub(":","_",substring(sequence_id,29))) %>%
                        select(label), new_tip , by= "label") %>%
    arrange(., genus_name)
}
  


# create new tax label for ease of plotting in heatmap:
new_tip$conservative_tax_new=NA
new_tip[is.na(new_tip[,"id"]),"conservative_tax_new"]=""

# fill in row 2+
for (i in c(2:(nrow(new_tip)))) {
  if (!is.na(new_tip[i,"id"])) {
    new_tip[i,"conservative_tax_new"] <-
      ifelse(new_tip[i,rank_label]!=new_tip[i-1,rank_label] | is.na(new_tip[i-1,rank_label]), 
             new_tip[i,rank_label], NA)
  }
}

# fill in row 1
new_tip[1, "conservative_tax_new"] <- new_tip[1,rank_label]

new_tip[is.na(new_tip$id),]$id <-ncol(ab)+1

new_tip$dot <- ifelse(is.na(new_tip$conservative_tax_new), NA, "·")





# outgroup's index as column name in heatmap
og <-as.character(ncol(ab)+1)

# save the image
pdf(paste("Fig.2_",tax_name,"_hmp.pdf", sep=""), width=8, height=8)


# heatmap, using rooted new tree
par(oma=c(2,1,0,1))
ht<-
  ab %>% as.data.frame %>% mutate(!!og := 0) %>% # include an empty column for outgroup
  select(new_tip$id)%>% as.matrix %>% # reorder matrix based on tree, since heatmap.2 doesn't do it automatically
  heatmap.2(hclustfun=function(x) hclust(x, method = meth),
            labRow = NA, 
            labCol = new_tip$dot, srtCol = 90, # labCol=new_tip$conservative_tax_new or $dot
            trace="none", key=T,
            col=hcl.colors(12, "Oslo", rev = F)
            ,Colv = NA,Rowv=NA, # ,Colv = drt # remove ,Rowv=NA to let rows cluster
            cexCol=0.8)# + geom_label_repel()

dev.off()

# remove intermediate files
system(paste("rm ",tax_name, "_new.fasta ", tax_name, "_new.nex", sep=""))




