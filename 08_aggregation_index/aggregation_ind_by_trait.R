
# Calculate nearest neighbour distance to evaluate distribution patterns of taxa

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
w <- owin(poly=list(x=c(9.03465, 9.18026,	9.19743,	9.23384,	9.25991,	9.33528,	9.32038,	9.37610,
                        9.41080,	9.47531,	9.42772,	9.48473,	9.61523,	9.6798,	9.75123,	9.71826,
                        9.62828,	9.59875,	9.53968,	9.48885,	9.41261,	9.34049,	9.06292),
                    y=c(47.80435,	47.72905,	47.64855,	47.63250,	47.64715,	47.62857,	47.59713,	47.57352,
                        47.61780,	47.5546,	47.52606,	47.48707,	47.53488,	47.50239,	47.50889,	47.55112,
                        47.56967,	47.58682,	47.59331,	47.64612,	47.66742,	47.66279,	47.8207)))


plot(w)


####################################################################
# Run "./10_add_trait/euk_add_trait_BTL.R" to add life mode to the EUK motus table
# To continue you need objects: dt, data
####################################################################


# a table of coordinates of each extract/core/site
coor <- left_join(dt$pcrs %>% mutate(sample_id=rownames(.)),
                  dt$samples %>% mutate(sample_id=rownames(.)), 
                  by=c("sample_id","sediment")) %>% 
  select(sediment,sample_id, long, lat, site)

# assign numbers to prepare for coord shift
c <- coor %>% 
  group_by(site) %>% mutate(i=as.numeric(factor(sediment))-1) %>% 
  group_by(sediment) %>% mutate(j=as.numeric(factor(sample_id))-1)

# shift coordinates: 20 cm between cores at a site; 2 cm between extracts from same core
coor_s <- c %>% mutate(long = long+i*3e-6+j*2e-7, lat = lat+i*3e-6+j*2e-7) 




tmp_agg <- data; tax_name="euk"; rank_name="phylum_name"

# subset data to drop MOTUs appearring in fewer than n sites
n=2
tmp_new <- subset_metabarlist(tmp_agg, "motus", colSums((tmp_agg$reads>0)*1) >= n)






x <- coor_s$long
y <- coor_s$lat


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






pdf("Fig.2d_clarkevans_violin_lifemode.pdf", height=5, width=7)

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






