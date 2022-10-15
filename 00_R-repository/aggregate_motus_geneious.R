
library(metabaR)
library(seqRFLP)

suppressPackageStartupMessages(library(dendextend)) # to manage clustering dendrogram
library(ape)
library(tidytree)
library(tibble)
library(phylogram)


aggregate_motus_geneious <- function(dt, tax_name){
  
  # branch length to agg ASVs
  v = 0.007 
  
  # choose tree outgroup based on tax
  if (tax_name == "Copepoda"){
    og <- c(
      "outgroup", # Ameridae
      "ttaggcggggtacggccgcctcgcgcggccgtcttctcggatcaatcgcggagcttagctctccgtcgagacccgc"
      )
    #v = 13.3 # a diff branch length for COP
    v = 40
    
  } else if (tax_name =="Fungi"){
    og <- c("outgroup", # Cercozoa
            "ggaaaacttaccaggtccagacatagtaaggattgacagattgaagatctttcttgattctatgggtggtggtgcatggccgttcttagttggtggagtgatttgtctggttaattccgttaacgaacgagacctcaacctgctaaataggtccacgaattcctcgggattcgtgtggccttcttagagggactatcggtgatttagccgacggaagtttgaggca"
      )
    
  } else if (tax_name =="Cyanobacteria"){
      og <- c(
        "outgroup", # Acidobacterium
        "ggaaaccctgacgaagcaacgccgcgtggaggatgaaggccttcgggtcgtaaactcctgtcgactgggaagaatgcacctgacctaatacgtcggcgtgttgactgtaccggtggaggaagccacggctaactctgtgccagcagccgcggtaatacagaggtggcaagcgttgttcggaattactgggcgtaaagggcgcgtaggcggcccgttaagtcccgtgtgaaagcccccggctcaaccggggaacggcgcgggaaactggcaggcttgagttcgggagagggaagcggaatttcgggtgtagcggtgaaatgcgcagatatccgaaggaacaccggtggcgaaggcggcttcctggaccgacactgacgctgaggcgcgaaagctaggggagcaaac"
      )
  } else {
    og <- c("outgroup", # Chaetonotus
            "ggaaaactcacccggccaggacaccgtaaggattgacagactgagagctctttcttgattcggtgggtggtggtgcatggccgttcttagttggtggagcgatttgtctggttaattccgataacgaacgagactctagcctgctaaatagacgggcaatcccattgtgggttgacccgatttgcttcttagagggacaagtggcgttcctagccacgcgaaattgagc"
    )
    }
  
  
  
  # output as fasta file
  dataframe2fas(dt$motus %>% mutate(id=substring(rownames(.),29)) %>% # keep tail of sequence id
                select(id, sequence) %>%
                rbind(og), 
              file=paste(tax_name, ".fasta", sep = ""))
  
  
  # Go to Geneious to align reads* (method: Clustal Omega)
  # and build phylogeny (Tamura-Nei, Neighbour-joining). 
  # Output tree in nexus format, replace characters with _, dt.nex (default)
  
  system(paste(
    "geneious -i ",
    tax_name,
    ".fasta -w ./08_heatmaps/test.geneiousWorkflow -o ",
    tax_name,
    ".nex", sep=""))
  
  

  
  nx <- read.nexus(file = paste(tax_name,".nex",sep = "")) # :- replaced to _
  df <- as_tibble(nx) # transform to data.frame class
  
  # Method 1 (conservative): see collapse_tips_1.R
  # Method 2: collapse a whole clade where longest branch is shorter than a given value.
  
  
  
  # Functions needed:
  trace_back <- function(tip, df) {
    # Inputs are 1) a tip node number, 2) the phylogeny in data.frame class.
    # Output is a list of all its parent nodes directly connected to this tip.
    out=c(tip, get_parent(tip, df))
    while (out[length(out)] != out[length(out)-1]) {
      out <- c(out, get_parent(out[length(out)], df))
    }
    return(unique(out))
  }
  
  get_parent <- function(node, df) {
    # get direct parent of a node in df (phylogenetic tree in data.frame class)
    return(df[df$node==node,]$parent)
  }
  
  
  
  
  # use v as subtree distance
  # filter out nodes with any branch longer than v
  watch <- df %>% filter(branch.length < v) %>% select(parent) %>% unique
  discard <-  df %>% filter(branch.length >= v) %>% select(parent) %>% unique
  
  keep <- watch # initialise table containing nodes to keep
  
  for (node in discard$parent) {
    keep <- keep %>% filter(!parent %in% trace_back(node, df))
  }
  
  keep_tb <- df %>% filter(parent %in% keep$parent)
  
  
  # assign collapsing groups
  pool <- unique(keep_tb$parent) # unique parents
  keep_tb$group <- match(keep_tb$parent, pool) # first assign group using the index number in pool
  
  # then update group assignment for clades containing > 2 tips.
  br <- keep_tb %>% filter(is.na(label)) # branches, i.e. not tips
  for (i in c(1:nrow(br))) {
    keep_tb[keep_tb$parent==br$node[i],]$group <- br$group[i]
  }
  
  
  
  
  
  # assign taxa info
  motus <- dt$motus %>% mutate(label=gsub(":","_",substring(rownames(.),29))) %>%
    select(label, SCIENTIFIC_NAME)
  
  keep_tb_tax <-left_join(keep_tb, motus, by="label") %>% filter(!is.na(label))
  
  
  
  
  # include only grouped nodes where taxa assignment is same for all seq
  to_clps <- keep_tb_tax %>% group_by(group) %>% 
    summarise(nTaxa=length(unique(SCIENTIFIC_NAME))) %>%
    filter(nTaxa==1)
  
  
  # create an updated group variable and save in a new table
  keep_tb_tax_final <- keep_tb_tax %>% mutate(group2 = ifelse(group %in% to_clps$group, group, NA)) %>% select(label, group2)
  
  # assign group numbers to all tips
  df2 <- left_join(df %>% filter(!is.na(label)), keep_tb_tax_final, by="label") %>%
    mutate(group3=ifelse(is.na(group2), row_number()+max(to_clps$group), group2))
  
  
  
  
  
  # create new motus table with assigned groups-to-collapse
  motus <- left_join(dt$motus %>% mutate(label=gsub(":","_",substring(rownames(.),29))),df2, by="label") %>% select(-c(label,group2))
  
  # final update group assignment so that group numbers are consecutive
  groups <- unique(motus$group3)
  motus$group_final <- match(motus$group3, groups); motus$group3<-NULL
  rownames(motus) <- motus$sequence_id
  
  # update motus table in the metabarlist
  dt$motus<- motus
  
  
  
  
  check_metabarlist(dt)
  dt2<- aggregate_motus(dt, groups=dt$motus$group_final)
 
  
  return(dt2)
}