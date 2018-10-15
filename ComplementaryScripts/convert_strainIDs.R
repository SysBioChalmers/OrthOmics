#Strains identifiers conversion: Gets a list of gene IDs for a given Y. lipolytica 
#strain (W29 or CLIB122) or K. marxianus (DMKU or CBS 6556) and returns a new one 
#with the other strain identifiers (limited those genes for which a single copy 
#ortholog is available)

convert_strainIDs <- function(org,strain,dataset){
  IDs <- dataset$genes
  if (all(tolower(org) == 'yli')){
    #Load annotation data: single copy orthologs W29 <-> CLIB122
    setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Orthologs')
    IDs_table <- read.delim('Yali_annotation.txt', header = TRUE, sep = "\t",stringsAsFactors=FALSE)
    #Discard internal CHASSY ids column
    IDs_table <- IDs_table[,2:3] 
    #Standardize format
    #IDs_table$CLIB122 <- gsub('YALI0','YALI0_',IDs_table$CLIB122)
  }
  
  if (all(tolower(org) == 'kma')){
    #Load annotation data: single copy orthologs DMKU <-> CBS6556
    setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Orthologs')
    IDs_table <- read.delim('Kma_annotation.txt', header = TRUE, sep = "\t",stringsAsFactors=FALSE)
    #Extract columns of interest
    IDs_table[,3] <- IDs_table[,6]
    IDs_table <- IDs_table[,2:3]
    colnames(IDs_table) <- c('DMKU','CBS')
  }
  
  #Identify strain column and get IDs
  column    <- match(tolower(strain), tolower(colnames(IDs_table)))
  strainIDs <- IDs_table[,column]
  newCol <- c(1,2)[c(1,2)!=column]
  #Get other strain's ids
  matches <- match(IDs,IDs_table[,column])
  #Get indexes for converted elements
  converted <- !is.na(matches)
  #Convert just those IDs for which a match was found
  newIDs  <- IDs
  newIDs[converted]  <- IDs_table[matches[converted],newCol]
  dataset$genes <- newIDs
  return(dataset)
}
