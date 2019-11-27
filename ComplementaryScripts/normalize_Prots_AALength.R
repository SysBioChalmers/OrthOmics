normalize_Prots_AALength <- function(data,genes,proteins,organism){
#normalize_Prots_AALength
#
#Function that gets a spectral counts proteomics dataset, divides the counts values
#for each protein by its AA chain length (queried from Uniprot). Then the transformed
#dataset is normalized by the sum of all values multiplied by their respective molecular
#weights. The final units are [mmol/g protein] for each entry.
#
# data      (dataframe) Spectral counts measurements for a set genes/proteins 
#           (rows) in n-conditions (columns).It should include all the replicates 
#           of the original dataset.
# genes     List of gene IDs for the dataset
# proteins  List of Uniprot protein IDs for the dataset
# organism  'kma' for K. marxianus; 'sce' for S. cerevisiae; 'yli' for Y. lipolytica
#
# output    list including: newData - Normalized dataset.The dataset filters out those proteins (rows), 
#                           for which either MW or AA chain length was not found in uniprot
#                           genes - list of genes in newData, filtering out those entries for which
#                           MW was not available.
#                           proteins - list of proteins in newData.
#
# Usage: outputLust <- normalize_Prots_AALength(data,genes,proteins,organism)
#
# Last modified: Ivan Domenzain. 2019-11-27
#
  
  filename <- paste('uniprot_',organism,'.txt',sep='')
  database <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  MWs      <- database[,6]
  SEQ      <- database[,7]
  if (all(organism != 'sce')){
    #For Kma and Yli the MW query is done based on DMKU and CLIB122
    #uniprot IDs respectively.
    DB_ids <- database[,1]
    ids    <- proteins
    if (all(organism == 'yli')){
      MWs <- database[,7]
      SEQ <- database[,8]
    }
  }else {
    #for sce the match is done based on CENPK gene IDs
    DB_ids <- database[,1]
    ids    <- genes
  }
  #Reduce dataset to those entries which have a correspondance in the
  #uniprot database
  DB_indxs <- match(ids,DB_ids)
  newData  <- c()
  newGenes <- c()
  newProts <- c()
  #Exclude indexes without a match in the database, with non-positive or
  #non-numerical MW
  MWs         <- as.numeric(MWs)/1000 #g/mmol
  avgMW       <- mean(MWs)
  NormFactors <- c()
  #get mean AA seq length
  AAlengths <- c()
  for (i in 1:length(DB_indxs)){
    index <- DB_indxs[i]
    if (!is.na(index)){AAlengths <- c(AAlengths,nchar(SEQ[index]))}
  }
  avgLength <- mean(AAlengths)
  
  for (i in 1:length(ids)){
    index <- DB_indxs[i]
    LengthAA <- avgLength
    Mweight  <- avgMW
    if (!is.na(index)){
      if ((!is.na(MWs[index]) | MWs[index]>0) & !is.na(SEQ[index])){
        LengthAA <- nchar(SEQ[index])
        Mweight  <- MWs[index]
      }
    }
    #Divide each protein row by its correspondant AA chain length
    newRow      <- data[i,]/LengthAA
    newData     <- rbind(newData,newRow)
    newRow      <- Mweight*newRow
    NormFactors <- rbind(NormFactors,newRow)
    newGenes    <- rbind(newGenes,genes[i])  
    newProts    <- rbind(newProts,proteins[i])  
  }
  NormFactors <- colSums(NormFactors)
  print(NormFactors)
  #Normalize dataset, brings the values back to their original order of
  #magnitude. Multiplying all columns by the same constant does not affect
  #fold-change calculations
  for (j in 1:ncol(newData)){
    newData[,j] <- newData[,j]/NormFactors[j]
  }
  
  return(list(newData,newGenes,newProts))
}