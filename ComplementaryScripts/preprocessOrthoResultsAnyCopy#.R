library(dplyr)
#library(stingr)

#preprocessOrthoResults <- function(repoPath){
#preprocessOrthoResultsAnyCopy#
#
#Function that generates a simple .csv with any orthologous proteins
#across "n" organisms from the output folder of orthoFinder
#  
# repoPath    (string) Directory path for this repository.
# newOGlist   (data frame) String containing each of the orthologous proteins in the test organism that match
#             any of the query organism(s).
#             
#
# Usage: newOGlist <- preprocessOrthoResults(repoPath)
#
# Last modified: Tyler Doughty 7-3-19
#
  
  #Open OrthoFinder results files
  #sourcePath <- paste(repoPath,'/OrthoFinder/Orthogroups',sep='')
  setwd('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Orthologs/orthfinder/s288cvgenus')
  filename      <- '/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Orthologs/orthfinder/s288cvgenus/Orthogroups/Orthogroups.GeneCount.tsv'
  anyCopyOG  <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  filename      <- '/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Orthologs/orthfinder/s288cvgenus/Orthogroups/Orthogroups.tsv'
  orthoGroups   <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  endCol        <- ncol(orthoGroups)
  newOGlist     <- c()
  group_number  <- 'GroupV.csv' #Add the group number for the analysis
  filename   <- '/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Orthologs/orthfinder/s288cvgenus/s288c_gene_prot.csv'
  gene_2_prot <- read.csv(filename, header = TRUE)


anyCopyOG <- filter(anyCopyOG, A_s288c !=0) #Remove orthogroups with 0 values for the query organism (these are orthologs between the other organisms but not the query)
#anyCopyOG  <- anyCopyOG [,c(1,3,2,4,5)]
anyCopyOG <- anyCopyOG[apply(anyCopyOG[c(3:5)],1,function(z) any(z!=0)),] #remove orthogroups that only have genes in the query organism

#Map each  OG to orthoGroups data frame
  for (i in 1:nrow(anyCopyOG)){
    OGid    <- anyCopyOG[i,1] 
    index   <- which(orthoGroups[,1]==OGid)
    rowData <- orthoGroups[index,2:endCol]
    newRow  <- OGid


    #Split  protein ID's for each organism
  for (j in 1:(endCol-1)){
      cellText  <- rowData[[j]]
      #cellText  <- strsplit(cellText,"")
      #cellText  <- cellText[[1]][2]
  
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
#newOGlist$newcolumn <- S+PNA
#newOGlist2 <- cbind(newOGsimple, c(1))
#newOGlist2 <- newOGlist2[,-1]

#newOGsplit <- c()
newOGsimple <- toString(newOGlist, width=NULL)
test <- as.list(strsplit(newOGsimple, ","))
num.el <- sapply(test, length)
newOGsimple <- cbind(unlist(test), rep(1:length(test), num.el))
colnames(newOGsimple) <- c("protein","orthogroup")
newOGsimple <- as.data.frame(newOGsimple)
newOGsimple$orthogroup <- gsub('1', 'S+EMP', newOGsimple$orthogroup)
newOGsimple$protein <- gsub('\\s+', '', newOGsimple$protein)

#Map each  protein to gene name
#for (i in 1:nrow(newOGsimple)){
m <- merge(newOGsimple, gene_2_prot, by.x = "protein", by.y = "protein", all.x= TRUE)
m  <- m [,c(1,3,2)]
#}
write.csv(m, file = group_number, row.names = TRUE)

#write.csv(newOGsimple, file = group_number, row.names = TRUE)
