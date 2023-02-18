
# Calculate nearest neighbour distance to evaluate distribution patterns of taxa
# simplified space
# Reference: Clark and Evans 1954

library(metabaR)
library(spatstat.core) # calculate clark and evans index, function: clarkevans
library(geosphere) # to use distm to calculate distances

library(dplyr)
library(data.table)
library(ggplot2)
library(colorspace)

source("./00_R-repository/aggregate_motus_geneious.R")

# poly window boundary for clark evans index, coord should be anticlockwise; check ?owin
w <- owin(poly=list(x=c(9,0,0,4,4,6,6,7,7,8,8,9),
                    y=c(5,5,4,4,0,0,1,1,0,0,3,3)))

plot(w)



# load data
dt <- readRDS("./00_R-repository/euk_r2_silva_97_agg_2rep_norm")
dt <- readRDS("./00_R-repository/cyano_r2_silva_97_agg_2rep_norm")




# add sequence_id into motus table, for later use
dt$motus <- dt$motus %>% mutate(sequence_id=rownames(.)) %>% select(sequence_id, everything())

# a table of coordinates of each extract/core/site
coor <- left_join(dt$pcrs %>% mutate(sample_id=rownames(.)),
                  dt$samples %>% mutate(sample_id=rownames(.)), 
                  by=c("sample_id","sediment")) %>% 
  select(sediment,sample_id, long, lat, site)

# assign idx to prepare for coord shift: same sediment with same i, diff extract diff j 
c <- coor %>% 
  group_by(site) %>% mutate(i=as.numeric(factor(sediment))-1) %>% 
  group_by(sediment) %>% mutate(j=as.numeric(factor(sample_id))-1) %>%
  mutate(i_old=i,j_old=j) %>% 
  mutate(i=ifelse(i_old==2&j_old==0,0,i), j=ifelse(i_old==2&j_old==0,1,j)) %>% # correct idx for Ueberlingen sites
  mutate(i=ifelse(i_old==3&j_old==0,1,i), j=ifelse(i_old==3&j_old==0,1,j)) %>%
  select(-i_old, -j_old)

 


# coordinates of sites
coor_st <- read.table("./11_author_response/site_coordinates_simplified_space.txt", header=T, sep='\t')

# shift coordinates:

coor_s <- c %>% left_join(coor_st, by="site") %>%
  mutate(extract_x = (site_x-0.25)-0.5*i, extract_y = (site_y-0.25)-0.5*j) 


# subset dataset
tax_name="Nematoda";rank_name = "phylum_name" 
tax_name="Chlorophyta";rank_name = "phylum_name" 
tax_name="Annelida";rank_name = "phylum_name" 
tax_name="Bacillariophyta";rank_name = "phylum_name" 
tax_name="Platyhelminthes";rank_name = "phylum_name" 
tax_name="Fungi";rank_name = "kingdom_name" 
tax_name="Arthropoda";rank_name = "phylum_name" 

tax_name="Cyanobacteria";rank_name = "phylum_name" 




tmp <- subset_metabarlist(dt, "motus",
                          dt$motus[,rank_name] == tax_name & !is.na(dt$motus[,rank_name]) &
                            !is.na(dt$motus$order_name))

# agg ASVs using phylogeny 
tmp_agg <- aggregate_motus_geneious(tmp, tax_name)

# subset data to drop MOTUs appearring in fewer than n sites
n=2
tmp_new <- subset_metabarlist(tmp_agg, "motus", colSums((tmp_agg$reads>0)*1)>=n)



x <- coor_s$extract_x
y <- coor_s$extract_y


if (!exists("out")){
  out <- data.frame(ind=NA, tax=NA, TAXID=NA)
}


for (i in c(1:ncol(tmp_new$reads))) {
  z <- (tmp_new$reads[,i]>0)*1
  X <- data.frame(x=x,y=y,z=z) %>% filter(z==1) %>% select(x,y) %>% unique
  pp <- as.ppp(X,w)
  
  out <- out %>% add_row(ind=clarkevans(pp, correction="cdf"),
                         tax=tax_name, TAXID=rownames(tmp_new$motus)[i]) 
}




tax_to_plot = c("Annelida" , "Platyhelminthes","Nematoda","Arthropoda","Bacillariophyta","Chlorophyta","Fungi", "Cyanobacteria")


pdf("ARFig.2c_clarkevans_violin_tax.pdf", height=5.5, width=8)

ggplot(out %>% filter(ind >= 0) %>% filter(tax %in% tax_to_plot), 
       aes(x=tax, y=ind, fill=tax)) +
  
  scale_x_discrete(limits = tax_to_plot) + 
  geom_violin(trim=T, scale="width", size=5, width=0.7, color=NA, adjust = 0.5, alpha=0.9) + 
  geom_jitter(shape=16, position=position_jitter(0.04), size=0.5, color="grey50") +
  geom_boxplot(width=0.2, size=0.4, alpha=0) + 
  
  scale_fill_brewer(palette="Set3") +
  ylab("Clark & Evans aggregation index") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1),
        axis.title.x=element_blank(), axis.ticks.x = element_blank(),
        legend.position="none",
        panel.grid = element_blank(),
        
        text = element_text(size=rel(4.6)),
        axis.title.y = element_text(margin = margin(r = 10)),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 

dev.off()



pdf("ARFig.2b_sampling_space_extracts.pdf", height=5, width=8)
plot(coor_s$extract_x, coor_s$extract_y, xlim=c(-1, 10), ylim=c(-1, 10), 
     xaxt='n',yaxt='n', bty = 'n', ann=F)
polygon(x=c(9,0,0,4,4,6,6,7,7,8,8,9),
        y=c(5,5,4,4,0,0,1,1,0,0,3,3))
dev.off()

pdf("ARFig.2b_sampling_space_sites.pdf", height=5, width=8)
plot(coor_s$extract_x, coor_s$extract_y, xlim=c(-1, 10), ylim=c(-1, 10)
     , col='white',xaxt='n', yaxt='n', bty = 'n', ann=F)
polygon(x=c(9,0,0,4,4,6,6,7,7,8,8,9),
        y=c(5,5,4,4,0,0,1,1,0,0,3,3))
clip(4,8,-100,100)
abline(h=c(1:3),lty='dashed')
clip(4,9,-100,100)
abline(h=4,lty='dashed')
clip(-100,100,4,5)
abline(v=c(1:4),lty='dashed')
clip(-100,100,0,5)
abline(v=c(5:8),lty='dashed')

dev.off()


