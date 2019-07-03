library(dplyr)
library(stingr)

preprocessOrthoResults <- function(repoPath){
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
  #setwd(sourcePath)
  filename      <- '/Users/doughty/Documents/GitHub/OrthoFinder-2.3.3_source/orthofinder/CBSvSCE-L.Ferm-E.goss/OrthoFinder/Results_Jun26/Orthogroups/Orthogroups.GeneCount.tsv'
  anyCopyOG  <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  filename      <- '/Users/doughty/Documents/GitHub/OrthoFinder-2.3.3_source/orthofinder/CBSvSCE-L.Ferm-E.goss/OrthoFinder/Results_Jun26/Orthogroups/Orthogroups.tsv'
  orthoGroups   <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  endCol        <- ncol(orthoGroups)
  newOGlist     <- c()
  group_number  <- 'GroupIII.csv' #Add the group number for the analysis
  

anyCopyOG <- filter(anyCopyOG, K.marxianus !=0) #Remove orthogroups with 0 values for the query organism (these are orthologs between the other organisms but not the query)
anyCopyOG  <- anyCopyOG [,c(1,3,2,4,5)]
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


newOGlist <- newOGlist[,-1:-2]
newOGlist <- newOGlist[,-2:-3]

newOGsplit <- c()
newOGsimple <- toString(newOGlist, width=NULL)
write.csv(newOGsimple, file = group_number, row.names = FALSE)
