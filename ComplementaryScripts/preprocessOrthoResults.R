preprocessOrthoResults <- function(repoPath){
#preprocessOrthoResults
#
#Function that generates a simple list with single copy orthologous proteins
#across "n" organisms from the output folder of orthoFinder
#  
# repoPath    (string) Directory path for this repository.
# newOGlist   (data frame) List containing each of the single copy orthologous
#             protein groups. OG IDs are placed in the first column and then uniprot
#             IDs for the proteins for each mapped organism in the subsequent columns.
#
# Usage: newOGlist <- preprocessOrthoResults(repoPath)
#
# Last modified: Ivan Domenzain. 2019-06-12
#
  
  #Open OrthoFinder results files
  sourcePath <- paste(repoPath,'/OrthoFinder/Orthogroups',sep='')
  setwd(sourcePath)
  filename      <- 'Orthogroups_SingleCopyOrthologues.txt'
  singleCopyOG  <- read.delim(filename, header = FALSE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  filename      <- 'Orthogroups.tsv'
  orthoGroups   <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  endCol        <- ncol(orthoGroups)
  newOGlist     <- c()
  #Map each single copy OG to orthoGroups data frame
  for (i in 1:nrow(singleCopyOG)){
    OGid    <- singleCopyOG[i,1]
    index   <- which(orthoGroups[,1]==OGid)
    rowData <- orthoGroups[index,2:endCol]
    newRow  <- OGid
    #Split  protein ID's for each organism
    for (j in 1:(endCol-1)){
      cellText  <- rowData[[j]]
      cellText  <- strsplit(cellText,"\\|")
      cellText  <- cellText[[1]][2]
      #Get a new row, binding together the OG id with the corresponding uniprot ID's 
      #for each of the queried organisms
      newRow    <- cbind(newRow,cellText)
    }
    newOGlist <- rbind(newOGlist,newRow)
  }
  colnames(newOGlist) <- colnames(orthoGroups)
  return(newOGlist)
}