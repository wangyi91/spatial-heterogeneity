# trait based heatmap for plants

library(metabaR)
library(dplyr)
library(data.table)
library(gplots) # to use heatmap.2
library(reshape2) # to use melt
library(ggplot2)

library(RColorBrewer) # to visualise color palette
library(dendextend) # to manage clustering dendrogram
library(colorspace)# to work with dendextend

source("./00_R-repository/aggregate_to_site.R")
source("./10_add_trait/plant_merge_traits_to_metabarlist.R")

# load data and transform to frequency of detection data at each site
dt_clean <- readRDS("./00_R-repository/plants_r2_embl_98_clean")
dt_toSite <- aggregate_to_site(dt_clean, stats="FoD")

# add traits to data and keep taxa that are identified to species and present in the trait table
dt <- plant_merge_traits_to_metabarlist(dt_toSite) %>%
  subset_metabarlist(., "motus", !is.na(.$motus$species_name) & !is.na(.$motus$info_flora))





ab <- dt$reads


# method for hclust
method="complete"

# other options
# method="ward.D"
# method="ward.D2"
# method="single"
# method="average" 
# method="mcquitty"
# method="median"
# method="centroid"





# ordinary heatmap
ht <- ab %>% 
  heatmap.2(hclustfun=function(x) hclust(x, method = method),
            trace="none", # labRow = NA, 
            col=c("#FCFCFC",hcl.colors(7, "Greens2", rev = T), "#004616"),
            Rowv=NA) # add Rowv=NA to disallow row clustering

# use the order in dendrogram for later use
ddorder <- ht$colInd

# for manually reorder dend to make a visually better gradient:
# make a column dendrogram just to see the size of each group
dd <- hclust(dist(t(ab)), method=method)

cutree(dd, k=3) %>% table
plot(dd, label=F)


# use the order in dendrogram for later use
# inspect and adjust the numbers each time, if needed
ddorder <- c(ht$colInd[248:303], ht$colInd[1:247])

# Arrange plant species differently for heatmap:
# 1. by aquatic, alpine, cultivated, others (Fig.3a)
# 2. by growth form (Supp.Fig.6c)

# choose which type of trait to plot
trait_type <- "traits"
trait_type <- "growth_form"


motus_sorted <- dt$motus %>% mutate(idx=seq.int(nrow(.))) %>% 
  slice(ddorder) %>% #  # first order rows as in dendrogram so that it is kept later
  select(idx,TAXID, family_name, genus_name, species_name, growth_form, traits) %>%
  arrange(get(trait_type))



# create new tax label for ease of plotting in heatmap:
motus_sorted$label=NA

# fill in row 2+
for (i in c(2:(nrow(motus_sorted)))) {
  motus_sorted[i,"label"] <-
    ifelse(motus_sorted[i,trait_type] != motus_sorted[i-1,trait_type], 
           motus_sorted[i,trait_type], NA)
  }

# fill in row 1
motus_sorted[1, "label"] <- motus_sorted[1,trait_type]
motus_sorted$dot <- ifelse(is.na(motus_sorted$label), NA, "·")



# set figure number in manuscript
fig_no = ifelse(trait_type == "traits", 'Fig.3a_plants_', 'Supp.Fig.6c_plants_')

# Plot
pdf(paste(fig_no, trait_type,'_hmp.pdf', sep=""), width=8, height=8)
par(oma=c(2,1,0,1))
ht <-
  ab %>% as.data.frame %>% select(motus_sorted$idx) %>% as.matrix %>%
  heatmap.2(hclustfun=function(x) hclust(x, method = method), dendrogram = "none",
            #lhei=lhei, lwid=lwid, lmat = lmat,
            trace="none", labCol = motus_sorted$dot, labRow = NA, 
            col=c("#FCFCFC",hcl.colors(9, "PuBu", rev = T)),
            #col=c("#FCFCFC", hcl.colors(7, "Greens2", rev = T), "#004616")
            Colv = NA,Rowv=NA) # add Rowv=NA to disallow row clustering
dev.off()



####### Supp.Fig.6b #######
for (trait in c("alpine","aquatic","cultivated","others")) {
  taxids <- motus_sorted$TAXID[motus_sorted$traits==trait]
  sub <- subset_metabarlist(dt, "motus", dt$motus$TAXID %in% taxids)
  det <- reshape2::melt((sub$reads>0)*1) %>% dplyr::rename(Site=Var1) %>% group_by(Site) %>% dplyr::summarise(ntaxa=sum(value))
  frq <- reshape2::melt(sub$reads) %>% filter(value>0); colnames(frq) <- c("Site","TAXID","fod")
  a <- ggplot() +
    geom_point(data=frq, aes(x=factor(Site, levels=rev(unique(Site))), y=fod),
               position="jitter", color="black", alpha=0.4) +
    ylim(c(0,1)) +
    ylab("Frequence of detection") +
    ggtitle(trait) +
    theme_classic() +
    theme(axis.line.x =element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          text = element_text(size=10))
    
  b <- ggplot() +
    geom_col(data=det,aes(x=factor(Site, levels=rev(unique(Site))), y=ntaxa)) +
    xlab("Site") + 
    ylab("Number of species detected") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          text = element_text(size=10))
  
  plt <- cowplot::plot_grid(a, b, align = "v", axis = "blr", ncol = 1, rel_heights = c(0.4, 0.6))
   
  if (trait=="alpine"){ab1 <- plt} 
  else if (trait=="aquatic"){ab2 <- plt} 
  else if (trait=="cultivated") {ab3 <- plt} 
  else {ab4 <- plt}
  
}

pdf("Supp.Fig.6b_plant_detection.pdf", width = 12, height = 8)
cowplot::plot_grid(ab1, ab2,ab3,ab4, labels = NULL, nrow=2)
dev.off()


