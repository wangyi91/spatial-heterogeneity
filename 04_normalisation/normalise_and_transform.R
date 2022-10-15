# normalisation of pcr read read count using DESeq2


library(metabaR)
library(DESeq2) # for normalisation
library(genefilter) # provides the shorth function
library(vegan) # decostand function
library(ggplot2)
library(dplyr)


# load data
dt_name <- "plants_r2_embl_98"
dt_name <- "euk_r2_silva_97"
dt_name <- "cyano_r2_silva_97"
dt_name <- "cop_r5_EmblSilvaCustom_85"

# choose format of data to load; choose transformation method
#dt_format <- "_clean"
dt_format <- "_agg"; tr <- "_lgchord"
dt_format <- "_agg_2rep"
dt_format <- "_agg_2rep_motuAgg"
dt_format <- "_agg_motuAgg"

dt<-readRDS(paste("./00_R-repository/",dt_name, dt_format, sep=""))


# 1. Normalisation of reads table with DESeq

# create DESeq data object
cds <- DESeqDataSetFromMatrix(countData = t(dt$reads), colData = dt$pcrs, ~1, tidy = FALSE)

# Scale the data by estimating size factors using the shorth function 
cds = estimateSizeFactors(cds, type="poscounts", locfunc = genefilter::shorth)

# Return the scaled counts from the previous step 
ncounts<-counts(cds, normalized=TRUE)

# create a new metabarlist for normalised dataset: dt_norm
dt_norm <- dt
dt_norm$reads <- t(ncounts)

rm(cds);rm(ncounts)


# save RDS
saveRDS(dt_norm, file = paste("./00_R-repository/", dt_name, dt_format, "_norm", sep=""))



# 2. Log or log-chord transformation:
  
# choose which transformation to do:
#tr <- "_lg"
#tr <- "_lgchord"

# compute log or log-chord transformed reads table and save it as new metabarlist "_tr"
dt_norm_tr <- dt_norm

if (tr == "_lg"){
  dt_norm_tr$reads <- log1p(dt_norm$reads) # log transform
} else {
  dt_norm_tr$reads <- decostand(log1p(dt_norm$reads), "norm") # log-chord transform (for PCA)
}

# save RDS
saveRDS(dt_norm_tr, file = paste("./00_R-repository/", dt_name, dt_format, "_norm", tr, sep=""))


