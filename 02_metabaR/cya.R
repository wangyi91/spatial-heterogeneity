# Data processing of cyano_silva_97

# Original tutorial of metabaR please check: 
# https://metabarfactory.github.io/metabaR/articles/metabaRF-vignette.html


library(metabaR)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)

## Preparation
# Load input tables as metabarlist into metabaR
dt <- tabfiles_to_metabarlist(file_pcrs = "./03_metabaR/cyano_pcrs.txt",
                              file_samples = "./03_metabaR/surfsedi_samples_corr.txt",
                              file_reads = "./03_metabaR/cyano_r2_silva_95_reads_modified.txt",
                              file_motus = "./03_metabaR/cyano_r2_silva_95_motus_modified.txt")



# Detection of contaminants from DNA extraction step:
dt <- contaslayer(dt, control_types = "extraction", output_col = "not_extr_conta")

# Detection of contaminants from PCR NTCs:
dt <- contaslayer(dt, control_types = "pcr", output_col = "not_pcr_conta")

# Compute relative abundance of all pcr contaminants together 
a <- data.frame(conta.relab = (rowSums(dt$reads[,!dt$motus$not_extr_conta]) + 
                                 rowSums(dt$reads[,!dt$motus$not_pcr_conta])) / 
                  rowSums(dt$reads))

# Add information on control types
a$control_type <- dt$pcrs$control_type[match(rownames(a), rownames(dt$pcrs))]


# flag pcrs with total contaminant relative abundance > 10% of reads)
dt$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(dt$pcrs), rownames(a))]>1e-1,  F, T)


# Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa 
dt$motus$target_taxon <- dt$motus$phylum_name=="Cyanobacteria" & !is.na(dt$motus$phylum_name)


# Add a column "not_degraded"
dt$motus$not_degraded <-
  ifelse(dt$motus$BEST_IDENTITY < 0.97, F, T)


# Flag pcrs with an acceptable (T) or unacceptable sequencing depth (F)
dt$pcrs$seqdepth_ok <-ifelse(rowSums(dt$reads) <= 0, F,T)


# Extract pcrs that has >0 read number and are not controls, for analysis
sub <- subset_metabarlist(dt, 
                          table="pcrs", 
                          indices = rowSums(dt$reads)>0 & 
                            dt$pcrs$type=="sample")


# Flag pcrs that are good relicates among each other
sub <- pcrslayer(sub, output_col = "replicating_pcr", plot=F)




# Visualise the pcr replicates and outliers behavour in PCA with function `check_pcr_repl`,
# and Distinguish between pcrs obtained from different samples
mds = check_pcr_repl(sub, groups=sub$pcrs$sample_id,
                     funcpcr = sub$pcrs$replicating_pcr)
mds + labs(color="extracts") + scale_color_discrete()


# report the flagging 
dt$pcrs$replicating_pcr <- NA
dt$pcrs[rownames(sub$pcrs),"replicating_pcr"] <- sub$pcrs$replicating_pcr


## Lowering tag-jumps
# Define a vector of thresholds to test
thresholds <- c(0,1e-3,2e-3, 1e-2,2e-2,3e-2,4e-2)

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(dt,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- dt$pcrs$control_type[match(tmp$sample, rownames(dt$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) + 
  geom_boxplot(color="grey40") + 
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.001"), col="orange", lty=2) + 
  geom_jitter(aes(color=controls), width = 0.15, alpha=0.3) + 
  scale_color_manual(values = c("brown", "cyan4"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=3) + 
  theme_bw() + 
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text.x = element_text(face="bold"),
        axis.text.x = element_text(angle=40, h=1), 
        legend.position = "none")



# Data cleaning and aggregation
# Use tag-jump corrected metabarlist with threshold identified above
tmp <- tests[["t_0.001"]]

# Filtering on motus: keep motus that are defined as TRUE in all the criterion below
tmp <- subset_metabarlist(tmp, "motus", 
                          indices = rowSums(tmp$motus[,c("target_taxon", "not_extr_conta","not_degraded")])==3)


# Subset on pcrs and exclude controls 
dt_clean <- subset_metabarlist(tmp,"pcrs", 
                               indices = tmp$pcrs$type == "sample" & 
                                 rowSums(tmp$pcrs[,c("low_contamination_level", 
                                                     "seqdepth_ok", "replicating_pcr")]) == 3 &
                                 is.na(tmp$pcrs[,c("low_contamination_level")])==F)


# Now check if previous subsetting leads to any empty pcrs or MOTUs
if(sum(colSums(dt_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(dt_clean$reads)==0)>0){
  print("empty pcrs present")
  dt_clean <- subset_metabarlist(dt_clean, "reads",rowSums(dt_clean$reads)>0)}


############ Visualisation: compare data before and after cleaning ############

# Get original data only for samples
tmp <- subset_metabarlist(dt, table = "pcrs",
                          indices = dt$pcrs$type == "sample" & rowSums(dt_clean$reads)>0)

# Add sediment sample id for checks
tmp$pcrs$sediment <- tmp$samples$sediment[match(tmp$pcrs$sample_id, rownames(tmp$samples))]

dt_clean$pcrs$sediment <-
  dt_clean$samples$sediment[match(dt_clean$pcrs$sample_id,
                                  rownames(dt_clean$samples))]

# Build PCoA ordinations 
mds1 <- check_pcr_repl(tmp, groups = tmp$pcrs$sediment)

mds2 <- check_pcr_repl(dt_clean, groups = dt_clean$pcrs$sediment)

# Custom colors
a <- mds1 + labs(color = "Sediment sample") + 
  scale_color_discrete() +
  theme(legend.position = "none") + 
  ggtitle("Raw data")

b <- mds2 + labs(color = "Sediment sample") +
  scale_color_discrete() + 
  ggtitle("Clean data")

# Assemble plots
leg <- get_legend(b + guides(shape=F) + 
                    theme(legend.position = "right", 
                          legend.direction = "vertical"))
ggdraw() +
  draw_plot(a, x=0, y=0, width = 0.38, height = 1) + 
  draw_plot(b + guides(color=F, shape=F), x=0.38, y=0, width = 0.38, height = 1) +
  draw_grob(leg, x=0.39, y=0)


####################################################################################





# remove BE124 and BE127 because they are redundant of S25
# remove BE038 and BE039 because they are redundant of S06
dt_clean <- subset_metabarlist(dt_clean, "samples", !rownames(dt_clean$samples) %in% c("BE124", "BE127","BE038", "BE039"))

# Save the cleaned and aggregated datasets for future reference
saveRDS(dt_clean, file = "./00_R-repository/cyano_r2_silva_97_clean")

