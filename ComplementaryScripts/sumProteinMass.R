sumProteinMass <- function(data,proteins,genes,organism){
#sumProteinMass
#
#Function that gets a semi-absolute proteomics dataset in and calculates the 
#total mass of protein, multiplying the abundance of each protein by its molecular
#weigth and getting a total sum at the end.
#
# data      (dataframe) Semi-absolute quantification data [mmol/ g Protein]
# proteins  List of uniprot IDs for the dataset
# organism  'kma' for K. marxianus; 'sce' for S. cerevisiae; 'yli' for Y. lipolytica
# o
#
# output    totalMass - Vector including the total mass of protein
#           [g measured protein / g of Total protein] for each sample in the dataset.
#
# Usage:  output <- sumProteinMass(data,IDs,organism)
#
# Last modified: Ivan Domenzain. 2019-04-26
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
      MWs <- database[,7]
    }
  }else {
    #for sce the match is done based on CENPK gene IDs
    DB_ids <- database[,3]
    ids    <- genes
  }
  #Reduce dataset to those entries which have a correspondance in the
  #uniprot database
  DB_indxs  <- match(ids,DB_ids)
  protsMass <- c()
  #Exclude indexes without a match in the database, with non-positive or
  #non-numerical MW
  MWs  <- as.numeric(MWs)/1000 #g/mmol
  for (i in 1:length(DB_indxs)){
    index <- DB_indxs[i]
    if (!is.na(index)){
      if (!is.na(MWs[index]) | MWs[index]>0){
        #Multiply each abundances row by its correspondant MW
        newRow    <- MWs[index]*data[i,]
        protsMass <- rbind(protsMass,newRow)
      }
    }
  }
  protsMass <- colSums(protsMass)
  print(protsMass)
  return(protsMass)
}