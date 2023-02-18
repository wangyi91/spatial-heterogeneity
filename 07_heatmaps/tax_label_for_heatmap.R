# Create label of heatmaps for Fig.2 and Fig.3A

library(ggrepel)
library(ggthemes)
library(ggplot2)
library(cowplot)



# If working on EUK, CYA and COP datasets, run this script after heatmap_euk_cya_cop.R
# Object needed: new_tip

exp <- new_tip %>% select(conservative_tax_new, dot)
exp$x <- 10
exp$y <- ifelse(is.na(exp$conservative_tax_new), NA, nrow(exp) - rownames(exp) %>% as.numeric())


pdf(paste("Fig.2_", tax_name,"_label.pdf", sep=""), width=10, height=10)

# set size of labels
if (tax_name=="Fungi") {
  force=7;size=6
} else if (tax_name=="Cyanobacteria") {
  force=7;size=7
} else if (tax_name=="Ciliophora") {
  force=4;size=7
} else {
  force=9;size=9
}

ggplot(exp, aes(x=x, y=y)) + geom_point(size=2, pch=95) + 
  geom_text_repel(data=exp, aes(label=conservative_tax_new),
                  force=force, point.padding=0.4, 
                  vjust=0.1,
                  direction='y',
                  size=size,
                  box.padding = 0.1,
                  nudge_x=-0.006, ylim=c(0,nrow(exp)),
                  segment.size=0.3) +
  theme_nothing()

dev.off()



# If working on PLANT dataset, run scripts below after running heatmap_plant_by_trait.R
# Object needed: motus_sorted, trait_type


exp <- motus_sorted %>% select(label, dot)
exp$x <- 10
exp$y <- ifelse(is.na(exp$label), NA, nrow(exp) - seq.int(nrow(exp)) %>% as.numeric())


pdf(paste("plants_", trait_type, "_label.pdf", sep=""), width=10, height=10)


force=9;size=20
ggplot(exp, aes(x=x, y=y)) + geom_point(size=2, pch=95) + 
  geom_text_repel(data=exp, aes(label=label),
                  force=force, point.padding=0.4, 
                  vjust=0.1,
                  direction='y',
                  size=size,
                  box.padding = 0.1,
                  nudge_x=-0.006, ylim=c(0,nrow(exp)),
                  segment.size=0.7) +
  theme_nothing()


dev.off()

