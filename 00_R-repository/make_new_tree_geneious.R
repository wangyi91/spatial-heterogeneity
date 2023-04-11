library(metabaR)
library(seqRFLP)


make_new_tree_geneious<- function(dt, tax_name) {
  # choose tree outgroup based on tax
  if (tax_name == "Copepoda"){
    og <- c(
      "outgroup", # Ameridae
      "ttaggcggggtacggccgcctcgcgcggccgtcttctcggatcaatcgcggagcttagctctccgtcgagacccgc"
    )
   
    
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
  dataframe2fas(dt$motus %>% mutate(id=gsub(":","_",substring(sequence_id,29))) %>% # keep tail of sequence id
                  select(id, sequence) %>%
                  rbind(og), 
                file=paste(tax_name, "_new.fasta", sep = ""))
  
  
  # Go to Geneious to align reads* (method: Clustal Omega)
  # and build phylogeny (Tamura-Nei, Neighbour-joining). 
  # Output tree in nexus format, replace characters with _, dt.nex (default)
  
  system(paste(
    "geneious -i ",
    tax_name,
    "_new.fasta -w ./07_heatmaps/run.geneiousWorkflow -o ",
    tax_name,
    "_new.nex", sep=""))
  
}
