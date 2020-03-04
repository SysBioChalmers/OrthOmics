#Function that generates a simple .csv with any orthologous proteins
#across "n" organisms from the output folder of orthoFinder
#FORCES DUPLICATE PROTEINS TO CO-MIGRATE INTO A SINGLE GROUP  
#Order - Genus (group IV), clade (group III), subphylum (group IV), phylum (group V)             

  setwd('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/GroupI')
  filename      <- 'Orthogroups.GeneCount.tsv'
  anyCopyOG  <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  filename      <- 'Orthogroups.tsv'
  orthoGroups   <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  endCol        <- ncol(orthoGroups)
  newOGlist     <- c()
  group_number  <- 'GroupIDUPEs_SCE.csv' 
  filename   <- '/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/s288c_gene_prot.csv'
  gene_2_prot <- read.csv(filename, header = TRUE)
  SelfSelfDuplicates <- read.delim('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/GroupV/Duplications_Groups_For_Sorting.csv', header = FALSE)
  
                    

multiCopyOG <- anyCopyOG[apply(anyCopyOG[c(2)],1,function(z) any(z!=0)),] #Organism of interest in column 2 - Remove orthogroups with 0 values for the query organism (these are orthologs between the other organisms but not the query)
multiCopyOG <- multiCopyOG[apply(multiCopyOG[c(2)],1,function(z) any(z>=2)),] # keep duplicated genes in SCE from this query

library(tidyverse)
multiCopyOG1 <- multiCopyOG %>% filter(multiCopyOG[c(2)]==multiCopyOG[c(3)])
multiCopyOG2 <- multiCopyOG %>% filter(multiCopyOG[c(2)]==multiCopyOG[c(4)])
multiCopyOG3 <- multiCopyOG %>% filter(multiCopyOG[c(2)]==multiCopyOG[c(5)])

merged <- rbind(multiCopyOG1, multiCopyOG2, multiCopyOG3, by=c(1))
merged <- merged[-nrow(merged),]
mergedunique <- unique(merged$Orthogroup)
mergedunique <- as.data.frame(mergedunique)
names(mergedunique)[names(mergedunique) == "mergedunique"] <- "Orthogroup"

#Map each  OG to orthoGroups data frame
  for (i in 1:nrow(mergedunique)){
    OGid    <- mergedunique[i,1] 
    index   <- which(orthoGroups[,1]==OGid)
    rowData <- orthoGroups[index,2:endCol]
    newRow  <- OGid


    #Split  protein ID's for each organism
  for (j in 1:(endCol-1)){
      cellText  <- rowData[[j]]
     
      #Get a new row, binding together the OG id with the corresponding uniprot ID's 
      #for each of the queried organisms
      newRow    <- cbind(newRow,cellText)
    }
    newOGlist <- rbind(newOGlist,newRow)
  }
  colnames(newOGlist) <- colnames(orthoGroups)
  return(newOGlist)


  #Get OG names for OGs with shared copy number in SCE and at least one other organism
newOGlist <- as.data.frame(newOGlist)
newOGlist <- newOGlist[c(2)]

#Reformat the Self-Self Duplicate list to show orthologous proteins as prot1, prot2, prot3, etc. This will be
#matched to the orthofinder result files
SelfSelfDuplicates <- gsub(", ,", "", SelfSelfDuplicates$V1, fixed = TRUE)
SelfSelfDuplicates <- gsub(" ,", "", SelfSelfDuplicates, fixed = TRUE)
SelfSelfDuplicates <- as.data.frame(SelfSelfDuplicates) 
SelfSelfDuplicates2 <- trimws(SelfSelfDuplicates$SelfSelfDuplicates, which = c("right"))
SelfSelfDuplicates <- as.data.frame(SelfSelfDuplicates2) 

colnames(SelfSelfDuplicates) <- "A_s288c"
library(dplyr)
#Create a list of Self-Self Duplicates that are not found in the search query
subtract1 <- (SelfSelfDuplicates %>% anti_join(newOGlist))
#Create a list of Self-Self Duplicates that are found in the search query
DupeOGlist <- (SelfSelfDuplicates %>% anti_join(subtract1))
  
#newOGsplit <- c()
newOGsimple <- toString(DupeOGlist$A_s288c, width=NULL)
test <- as.list(strsplit(newOGsimple, ","))
num.el <- sapply(test, length)
newOGsimple <- cbind(unlist(test), rep(1:length(test), num.el))
colnames(newOGsimple) <- c("protein","Group")
newOGsimple <- as.data.frame(newOGsimple)
newOGsimple$Group <- gsub('1', 'Phylum-Duplicates-SCE', newOGsimple$Group)
newOGsimple$protein <- gsub('\\s+', '', newOGsimple$protein)

#Map each  protein to gene name

m <- merge(newOGsimple, gene_2_prot, by.x = "protein", by.y = "protein", all.x= TRUE)
m  <- m [,c(1,3,2)]
setwd('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting')
write.csv(m, file = group_number, row.names = FALSE)
