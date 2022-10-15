# trait based heatmap for *plants*.

library(vegan)
library(metabaR)
library(dplyr)
library(data.table)
library(gplots) # to use heatmap.2

library(RColorBrewer) # to visualise color palette
library(dendextend) # to manage clustering dendrogram
library(colorspace)# to work with dendextend

# on abundance data
dt <- readRDS("./00_R-repository/plants_r2_embl_98_agg_2rep_motuAgg_norm_sp_trait"); ab <- log1p(dt$reads)
# or freq of detection data - Use this!!
dt <- readRDS("./00_R-repository/plants_r2_embl_98_motuAgg_FoD_trait")
dt <- subset_metabarlist(dt, "motus", !is.na(dt$motus$species_name) & !is.na(dt$motus$info_flora)); ab <- dt$reads


# method for hclust
method="ward.D2"

# less good options
method="ward.D"
method="complete"
method="single"
method="average" 
method="mcquitty"
method="median"
method="centroid"





# ordinary heatmap
ht<- ab %>% 
  heatmap.2(hclustfun=function(x) hclust(x, method = method),
            trace="none", # labRow = NA, 
            col=c("#FCFCFC",hcl.colors(7, "Greens2", rev = T), "#004616"),
            Rowv=NA) # add Rowv=NA to disallow row clustering

# use the order in dendrogram for later use
ddorder <- ht$colInd

# for manually reorder dend to make a visually better gradient:
# make a column dendrogram just to see the size of each group
dd <- hclust(dist(t(ab)), method=method)

cutree(dd, k=2) %>% table
plot(dd, label=F)


# use the order in dendrogram for later use
# inspect and adjust the numbers each time
ddorder <- c(ht$colInd[101:131], ht$colInd[1:100])














# Arrange plant  species based a few ways:
# 1. by growth form
# 2. by native, cultivated, alpine, aquatic


# choose which type of trait to plot
trait_type <- "growth_form"
trait_type <- "traits"

# 
motus_sorted <- dt$motus %>% mutate(idx=seq.int(nrow(.))) %>% 
  dplyr::slice(all_of(ddorder)) %>% #  # first order rows as in dendrogram so that it is kept later
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




# save plot to local
pdf(paste('Fig.3A_plants_',trait_type,'_hmp.pdf', sep=""), width=8, height=8)


# heatmap
par(oma=c(2,1,0,1))

# paste code here
ht<-
  ab %>% as.data.frame %>% select(motus_sorted$idx) %>% as.matrix %>%
  heatmap.2(hclustfun=function(x) hclust(x, method = method), dendrogram = "none",
            #lhei=lhei, lwid=lwid, lmat = lmat,
            trace="none", labCol = motus_sorted$dot, labRow = NA, 
            col=c("#FCFCFC",hcl.colors(9, "PuBu", rev = T)),
            #col=c("#FCFCFC", hcl.colors(7, "Greens2", rev = T), "#004616")
            Colv = NA,Rowv=NA) # add Rowv=NA to disallow row clustering
dev.off()
































**Manually examine and modify each time!**
  ```{r}
# use the order in dendrogram for later use
# motuAgg:
ddorder <- c(ht$colInd[264:829], ht$colInd[1:263])
```




Categorize grown forms into groups
```{r}
data$motus <- data$motus %>% mutate(growth_form = ifelse(PlantGrowthForm %like% paste(c("tree" ,"shrub"), collapse="|"),
                                                         "tree/shrub",
                                                         PlantGrowthForm)) 
```



```{r}
motus_sorted <- data$motus %>% mutate(id=c(1:nrow(.))) %>% # give each ASV an id number
  slice(all_of(ddorder))# %>% # first order rows as in dendrogram so that it is kept later
select(id, growth_form, PlantGrowthForm, Aquatic, class_name, order_name, family_name, genus_name, species_name, SCIENTIFIC_NAME) %>%
  arrange(., growth_form)#, order_name, family_name, genus_name, species_name, SCIENTIFIC_NAME)


ht<-
  ab %>% as.data.frame %>% 
  select(motus_sorted$id)%>% as.matrix %>% # reorder matrix based on tree, since heatmap.2 doesn't do it automatically
  heatmap.2(hclustfun=function(x) hclust(x, method = 'ward.D2'),
            labRow = NA, 
            labCol = motus_sorted$PlantGrowthForm, srtCol = 90, # labCol=new_tip$conservative_tax_new or $dot
            trace="none",key=T,
            col=c("#FCFCFC",hcl.colors(7, "Greens2", rev = T), "#004616"),
            Colv = NA,Rowv=NA, # ,Colv = drt # remove ,Rowv=NA to let rows cluster
            cexCol=0.8)# + geom_label_repel()

```














```{r}
# get the matrix of heatmap and add identification
mat <- ht$carpet %>% as.data.frame %>% mutate(TAXID=as.numeric(rownames(.))) %>% left_join(., data$motus %>% select(TAXID, SCIENTIFIC_NAME), by="TAXID")
```

```{r}
# print names of selected rows of taxa
start <- which(grepl(3319, mat$TAXID))
end <- which(grepl(4022, mat$TAXID))

names <- mat %>% select(TAXID, SCIENTIFIC_NAME) %>% slice(start:end) %>% select(SCIENTIFIC_NAME)

write.table(names, "/Users/yiwang/Downloads/names.txt", sep="\n",row.names =F,col.names =F, quote = F)
```



```{r}
# prepare new tree tip labels: change nx_new_rt$tip.label to colnames of "ab", so that heatmap can use the tree in plotting
new_tip <- left_join(tree$tip.label %>% as.data.frame %>% `colnames<-`("SCIENTIFIC_NAME"), 
                     data$motus %>% mutate(id=rownames(.)) %>% select(SCIENTIFIC_NAME,id),
                     by="SCIENTIFIC_NAME")

nx_new_rt$tip.label<-new_tip$id



# dendrogram of taxonomy tree
dtx <-phylogram::as.dendrogram.phylo(tree)
labels(dtx)<-NULL;plot(dtx)
```



```{r}
# heatmap, using phylo tree
ht<-
  ab %>% as.data.frame %>% 
  select(new_tip$id)%>% as.matrix %>% # reorder matrix based on tree, since heatmap.2 doesn't do it (no idea why)
  heatmap.2(hclustfun=function(x) hclust(x, method = 'average'),
            labRow = NA, labCol = NA, trace="none",
            col=hcl.colors(12, "Oslo", rev = F)
            ,Colv = dtx,Rowv=NA) # ,Colv = drt # add ,Rowv=NA to disallow rows cluster
```























