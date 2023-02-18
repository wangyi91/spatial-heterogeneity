
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


####################################################################
# Run "./10_add_trait/euk_add_trait_BTL.R" to add life mode to the EUK motus table
# objects needed: dt, data
####################################################################

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


tmp_agg <- data; tax_name="euk"; rank_name="phylum_name"

# subset data to drop MOTUs appearring in fewer than n sites
n=2
tmp_new <- subset_metabarlist(tmp_agg, "motus", colSums((tmp_agg$reads>0)*1)>=n)



x <- coor_s$extract_x
y <- coor_s$extract_y


if (!exists("out")){
  out <- data.frame(ind=NA, tax=NA)
}


for (i in c(1:ncol(tmp_new$reads))) {
  z <- (tmp_new$reads[,i]>0)*1
  X <- data.frame(x=x,y=y,z=z) %>% filter(z==1) %>% select(x,y) %>% unique
  pp <- as.ppp(X,w)
  
  out <- out %>% add_row(ind=clarkevans(pp, correction = "cdf"),
                         tax=tmp_new$motus$Lebensform[i]) 
  
}





pdf("ARFig.2d_clarkevans_violin_trait.pdf", height=5.5, width=8)

ggplot(out %>% filter(ind >= 0), aes(x=tax, y=ind, fill=tax)) + #x=tax x=id
  
  scale_x_discrete(limits=c("P", "PB","B", "Pa", "E",  "T"),
                   labels=c("planktonic","planktonic / benthic",  "benthic", "parasitic","epiphytic / epizoic",  "terrestrial")) + #,"E"  "T"  "NS" "W"  "Em" "S"  "SE"
  
  geom_violin(trim=T, scale="area", size=5, width=0.7, color=NA, adjust = 0.5, alpha=0.9) + 
  geom_jitter(shape=16, position=position_jitter(0.04), size=0.5, color="grey50") +
  geom_boxplot(width=0.2, size=0.4, alpha=0) + 
  
  scale_fill_brewer(palette="Set3") +
  xlab("mode of life") + 
  ylab("Clark & Evans aggregation index") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1),
        axis.title.x=element_blank(), axis.ticks.x = element_blank(),
        axis.title.y=element_blank(), 
        legend.position="none",
        text = element_text(size=rel(4.6)),
        #axis.title.y = element_text(margin = margin(r = 10)),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 

dev.off()


