normalize_SCounts <- function(data,genes,proteins,organism){
#normalize_SCounts
#
#Function that gets a spectral counts proteomics dataset and normalizes it
#according to the individual protein MWeights (queried) from uniprot database
#files.
#
# data      (dataframe) Spectral counts measurements for a set genes/proteins 
#           (rows) in n-conditions (columns).It should include all the replicates 
#           of the original dataset.
# organism  'kma' for K. marxianus; 'sce' for S. cerevisiae; 'yli' for Y. lipolytica
#
# output    list including: newData - Normalized dataset with respect to the individual 
#                           MWs of the proteins.The dataset filters out those proteins (rows), 
#                           for which a MW was not found in uniprot
#                           genes - list of genes in newData, filtering out those entries for which
#                           MW was not available.
#                           proteins - list of proteins in newData.
#
# Usage: outputLust <- normalize_SCounts(data,genes,proteins,organism)
#
# Last modified: Ivan Domenzain. 2019-02-25
#
  
  filename <- paste('uniprot_',organism,'.txt',sep='')
  database <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  MWs      <- database[,6]
  if (all(organism != 'sce')){
    #For Kma and Yli the MW query is done based on DMKU and CLIB122
    #uniprot IDs respectively.
    DB_ids <- database[,1]
    ids    <- proteins
    if (all(organism == 'yli')){
      MWs    <- database[,7]
    }
  }else {
    #for sce the match is done based on CENPK gene IDs
    DB_ids <- database[,3]
    ids    <- genes
  }
  #Reduce dataset to those entries which have a correspondance in the
  #uniprot database
  DB_indxs   <- match(ids,DB_ids)
  toRemove   <- c()
  newData    <- c()
  #Exclude indexes without a match in the database, with non-positive or
  #non-numerical MW
  MWs <- as.numeric(MWs)
  for (i in 1:length(DB_indxs)){
    index <- DB_indxs[i]
    if (!is.na(index)){
      if (is.na(MWs[index]) | MWs[index]<=0){
        toRemove <- c(toRemove,i)
      } else{
        #Multiply each protein row by its correspondant MW
        newData <- rbind(newData,data[i,]*MWs[index])  
      }
    } else{
      toRemove <- c(toRemove,i)
    }
  }
  
  #Normalize dataset
  for (j in 1:ncol(newData)){
    newData[,j] <- newData[,j]/sum(newData[,j]) 
  }
  if (length(toRemove>1)){
    genes    <- genes[-toRemove]
    proteins <- proteins[-toRemove] 
  }
  return(list(newData,genes,proteins))
}