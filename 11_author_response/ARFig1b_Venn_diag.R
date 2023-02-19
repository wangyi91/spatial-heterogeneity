library(ggvenn)
library(metabaR)
library(cowplot)
library(dplyr)
#library(RColorBrewer) # to visualise color palette
library(ggplot2)
library(ggthemes)


# import data
#dt_name ="plants_r2_embl_98"
#dt_name ="euk_r2_silva_97"

for (dt_name in c("plants_r2_embl_98","euk_r2_silva_97","cyano_r2_silva_97","cop_r5_EmblSilvaCustom_85")) {
  dt <- readRDS(paste("./00_R-repository/", dt_name, "_agg_norm_lgchord", sep="")) %>%
    subset_metabarlist(., "samples", !.$samples$site %in% c("S01","S02","S06","S07","S09","S10","S11","S21","S22","S23","S24","S25"))
  
  # group east: "S03","S04","S05","S08"
  # group deepest: S14
  # group north&south: S12,S13,S15,S16
  # group west: S17,S18,S19,S20
  
  dt$pcrs <- dt$pcrs %>% mutate(group = ifelse(dt$samples[dt$samples$sediment==sediment,]$site %in% c("S14"), "deepest",
                   ifelse(dt$samples[dt$samples$sediment==sediment,]$site %in% c("S03","S04","S05","S08"), "east",
                          ifelse(dt$samples[dt$samples$sediment==sediment,]$site %in% c("S12","S13","S15","S16"), "north_south",
                                 "west"))))
  # group dt
  dt_gr <- dt
  dt_gr$pcrs$sample_id = dt_gr$pcrs$group
  dt_gr$samples <- data.frame(matrix(ncol=1,nrow=4))
  
  rownames(dt_gr$samples) <- c("north_south","west","deepest","east")
  
  colnames(dt_gr$samples) = "sites"
  dt_gr$samples$ghost <- NA
  
  gr <- aggregate_pcrs(dt_gr, replicates=dt_gr$pcrs$group, FUN_agg_pcrs_mean)
  tab <- gr$reads %>% t %>% as.data.frame

  
  # transform to list as input data of ggvenn
  x <- list(
    deepest = rownames(tab %>% filter(deepest>0)), 
    east = rownames(tab %>% filter(east>0)), 
    north_south = rownames(tab %>% filter(north_south>0)), 
    west = rownames(tab %>% filter(west>0))
  )
  
  v <- ggvenn(
    x, 
    fill_color = c("grey", "#CD534CFF", "olivedrab", "#0073C2FF"),
    stroke_size = 0.5, set_name_size = 3, text_size = 2,
    ) + theme(plot.title = element_text(hjust = 0.5, size = 13))

  
  if (startsWith(dt_name,"p")) {
    a1 <- v + ggtitle("Plant") 
  } else if (startsWith(dt_name,"e")) {
    a2 <- v + ggtitle("Eukatyote") 
  } else if (startsWith(dt_name,"cy")) {
    a3 <- v + ggtitle("Cyanobacteria") 
  } else {
    a4 <- v + ggtitle("Copepod") 
  }
  
  
}



# save figure
pdf("ARFig.1b_venn.pdf", width=8, height=6)
plot_grid(a1, a2, a3, a4, labels = NULL, nrow=2, scale = 1.1)
dev.off()

png("ARFig.1b_venn.png", width=8, height=6, units = 'in', res = 600)
plot_grid(a1, a2, a3, a4, labels = NULL, nrow=2)
dev.off()

