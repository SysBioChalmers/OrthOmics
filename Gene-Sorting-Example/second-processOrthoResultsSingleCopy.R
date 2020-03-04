
#Function that generates a simple .csv for single copy orthologous proteins
#across "n" organisms from the output folder of orthoFinder
#  

  setwd('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/GroupI')
  filename      <- 'Orthogroups.GeneCount.tsv'
  anyCopyOG  <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  filename      <- 'Orthogroups.tsv'
  orthoGroups   <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  endCol        <- ncol(orthoGroups)
  newOGlist     <- c()
  group_number  <- 'GroupI-Singles.csv' #Add the group number for the analysis
  filename   <- '/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/s288c_gene_prot.csv'
  gene_2_prot <- read.csv(filename, header = TRUE)
  DuplicatesToExclude <- read.csv('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/GroupV/Duplicates_to_Exclude.csv', header = TRUE,  sep = ",",stringsAsFactors=FALSE, na.strings = "NA")


singleCopyOG <- anyCopyOG[apply(anyCopyOG[c(2)],1,function(z) any(z!=0)),] #Organism of interest in column 2 - Remove orthogroups with 0 values for the organism of interest (these are orthologs between the other organisms but not the query)
singleCopyOG <- singleCopyOG[apply(singleCopyOG[c(2)],1,function(z) any(z<2)),] # remove duplicates for the organism of interest
singleCopyOG <- singleCopyOG[apply(singleCopyOG[c(3:5)],1,function(z) any(z!=0)),] #remove orthogroups that only have genes in the query organism


#Map each  OG to orthoGroups data frame
  for (i in 1:nrow(singleCopyOG)){
    OGid    <- singleCopyOG[i,1] 
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



newOGlist <- newOGlist[,-1]
newOGlist <- newOGlist[,-2:-4]


newOGsimple <- toString(newOGlist, width=NULL)
test <- as.list(strsplit(newOGsimple, ","))
num.el <- sapply(test, length)
newOGsimple <- cbind(unlist(test), rep(1:length(test), num.el))
colnames(newOGsimple) <- c("protein","orthogroup")
newOGsimple <- as.data.frame(newOGsimple)
newOGsimple$orthogroup <- gsub('1', 'S+PNA-Single', newOGsimple$orthogroup)
newOGsimple$protein <- gsub('\\s+', '', newOGsimple$protein)

#Remove proteins that were previously found to be duplicates (inside the DuplicatesToExclude list)
library(dplyr)
minusDuplicates <- (newOGsimple %>% anti_join(DuplicatesToExclude))

#Map each  protein to gene name
m <- merge(minusDuplicates, gene_2_prot, by.x = "protein", by.y = "protein", all.x= TRUE)
m  <- m [,c(1,3,2)]

setwd('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting')
write.csv(m, file = group_number, row.names = FALSE)
