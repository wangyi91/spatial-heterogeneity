library(metabaR)
library(dplyr)
library(reshape2)
library(stringr) # for string manipulation
library(gplots) # to use heatmap.2
library(RColorBrewer) # to visualise color palette
library(dendextend) # to manage clustering dendrogram
library(colorspace)# to work with dendextend


dt<-readRDS("./00_R-repository/plants_r2_embl_98_agg_norm")

# add sequence_id into motus table, for later use
dt$motus <- dt$motus %>% mutate(sequence_id=rownames(.)) %>% select(sequence_id, everything())

# aggregate motus
dt <- aggregate_motus(dt, groups=dt$motus$TAXID, FUN = FUN_agg_motus_sum)

# exclude bad species assignments
sub <- subset_metabarlist(dt, "motus", !dt$motus$SCIENTIFIC_NAME %like% paste(c("sp." ,"subgen.", "var","sect","subsect","Group") , collapse="|"))

# exclude taxa without family info:
sub <- subset_metabarlist(sub, "motus", (!is.na(sub$motus$family_name) & sub$motus$family_name=="Asteraceae") |
                            !is.na(sub$motus$genus_name))


# load data from info flora
IF <- read.csv("./10_add_trait/info_flora_Checklist_simplified.csv",header=T, na.strings="") %>% filter(Rangstufe == "sp")

# mark what species are present in info flora checklist 
sub$motus <- sub$motus %>% mutate(info_flora = sub$motus$species_name %in% paste(IF$Gattung, IF$Spezies, sep=" "))

# get a species list
IF_sp <- sub$motus[sub$motus$info_flora,]$species_name %>% unique

# non-cultivated plants: A, I, ni
ncul <- IF %>% filter(paste(IF$Gattung, IF$Spezies, sep=" ") %in% IF_sp) %>% filter(Indigenat.CH %in% c("I","A","I/N","A/N","ni"))

# cultivated plants: AC, NC
neo_cul <- IF %>% filter(paste(IF$Gattung, IF$Spezies, sep=" ") %in% IF_sp) %>% filter(Indigenat.CH %in% c(NA,"NC","AC","N"))


# subset data
# remove cultivated taxa
sub <- subset_metabarlist(sub,"motus", !sub$motus$species_name %in% paste(neo_cul$Gattung, neo_cul$Spezies,sep=" "))

# keep taxa registered in info flora
sub <- subset_metabarlist(sub,"motus", sub$motus$info_flora==T)


# clean up metabarlist
pcrs_keep <- rowSums(sub$reads) %>% as.data.frame %>% filter(.>0) %>% rownames
sub <- subset_metabarlist(sub, "pcrs", rownames(sub$pcrs) %in% pcrs_keep)
check_metabarlist(sub)


# Aggregate pcrs by site, output *detection probability* or *mean read count*.
dt_st <- sub
dt_st$pcrs <- left_join(sub$pcrs, sub$samples %>% mutate(sample_id=rownames(.)), by="sample_id")
rownames(dt_st$pcrs) <- rownames(sub$pcrs)
dt_st$pcrs$sample_id <- dt_st$pcrs$site

# load the new 'samples' table: sites
sites<- read.table("./02_metabaR/surfsedi_sites.txt", header=T, sep='\t')
rownames(sites) <- sites$site
dt_st$samples <- sites 


# get log mean read count at each site
data <- aggregate_pcrs(dt_st, FUN=FUN_agg_pcrs_mean)
data$reads <- log1p(data$reads) # log transform


# remove motus appearing in only 1 site, for better focused heatmap
motus_fewer2 <- colSums((data$reads>0)*1) %>% as.data.frame %>% filter(.<2) %>% rownames
data <- subset_metabarlist(data, "motus", !data$motus$TAXID %in% motus_fewer2)


# clustering method for hclust
method="complete"

#pa <- (data$reads>0)*1 # presence-absence 
ab <- data$reads # log read count

# for manually reorder dend to make a visually better gradient:
# make a column dendrogram just to see the size of each group
#dd <- hclust(dist(t(ab)), method=method)

#cutree(dd, k=2) %>% table
#plot(dd, label=F)



ht <-
  ab %>% 
  heatmap.2(hclustfun=function(x) hclust(x, method = method),
            trace="none", # labRow = NA, 
            col=c("#FCFCFC", hcl.colors(7, "Greens2", rev = T), "#004616"),
            Rowv=NA) # add Rowv=NA to disallow row clustering

# mark alpine taxa
#alp <- ifelse(data$motus[ht$colInd,]$alpine,data$motus[ht$colInd,]$SCIENTIFIC_NAME,NA)
# mark info flora non-cultivated flora
#ifncul <- ifelse(data$motus[ht$colInd,]$species_name %in% paste(ncul$Gattung, ncul$Spezies, sep=" "), data$motus[ht$colInd,]$SCIENTIFIC_NAME, NA) # data$motus[ht$colInd,]$SCIENTIFIC_NAME "--"

# mark info flora cultivated flora
#ifneo_cul <- ifelse(data$motus[ht$colInd,]$species_name %in% paste(neo_cul$Gattung, neo_cul$Spezies, sep=" "), "--", NA) #

# taxon names to label the heatmap
taxa <- ifelse(!is.na(data$motus$species_name), data$motus$species_name, data$motus$genus_name)



pdf('Supp.Fig.6a_plants_hmp_keySp.pdf', width = 13.5, height = 9)

lmat <- rbind(c(1,3,4), c(2,1,4))
lhei <- c(0.5, 1)
lwid <- c(1.5, 4, 0.75)
par(oma=c(5,4,4,2))

ht <-
  ab %>% 
  heatmap.2(hclustfun=function(x) hclust(x, method = method), dendrogram = "none",
            lhei=lhei, lwid=lwid,lmat = lmat,
            trace="none", labCol=as.expression(lapply(taxa, function(a) bquote(italic(.(a))))), # labRow = NA, 
            col=c("#FCFCFC", hcl.colors(14, "PuBuGn", rev = T)[2:14]),
            Rowv=NA) # add Rowv=NA to disallow row clustering

dev.off()









