load_ProtData <- function(dataPath,organism,Pmethod){
  setwd(dataPath)
  pointer <- which(c('kma','sce','yli') == organism)
  if (all(Pmethod == 'XIC')){
    filename       <- paste(organism,'_XIC.txt',sep='')
    scalingFactors <- c(1,1,1)
  }
  if (all(Pmethod == 'SCounts')){
    filename       <- paste(organism,'_SCounts.txt',sep='')
    scalingFactors <- c(1,1,1)
  }
  if (all(Pmethod == 'NSAF')){
    filename       <- paste(organism,'_NSAF_bulk_fmol.txt',sep='')
    scalingFactors <- c(800,800,983) #ng of peptides in the sample
  }
  if (all(Pmethod == 'IBAQ')){
    filename <- paste(organism,'_IBAQ_bulk_fmol.txt',sep='')
    scalingFactors <- c(492,750,470) #ng of peptides in the sample
  }
  dataset <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  # Add the grouping information (and enforce the order of the samples).
  if (all(organism == 'yli')){
    conditions <- c('Ref','HiT','LpH')
    colorValues <- c("black", "red", "#009E73")
    replicates <- c(4,3,4)
    group <- factor(c(rep(conditions[1],replicates[1]),
                      rep(conditions[2],replicates[2]),
                      rep(conditions[3],replicates[3])), levels = conditions)
  } else {
    conditions <- c('Ref','HiT','LpH','Osm')
    colorValues <- c("black", "red", "#009E73","blue")
    replicates <- c(3,3,3,3)
    group <- factor(c(rep(conditions[1],replicates[1]),
                      rep(conditions[2],replicates[2]),
                      rep(conditions[3],replicates[3]),
                      rep(conditions[4],replicates[4])), levels = conditions)
    }
  IDs      <- dataset[,1:3]
  dataset  <- dataset[,(ncol(IDs)+1):(sum(replicates)+ncol(IDs))]
  dataset[is.na(dataset)] <- 0
  #Normalize data fmol -> [fmol / ng protein]
  dataset  <- dataset/scalingFactors[pointer] 
  genes    <- IDs[,1]
  if (all(organism == 'sce')){
    genes    <- IDs[,3]
    genes[genes==''] <- IDs[,1][genes=='']
  }
  proteins <- IDs[,2]
  proteins[proteins==''] <- genes[proteins=='']
  #Generate colNames
  cols <- c()
  for (i in 1:length(conditions)) {
    for (j in 1:replicates[i]) {
      cols <- c(cols,paste(conditions[i],'_',j,sep=''))
    }
  }
  colnames(dataset) <- cols
  remove(cols)
  return(list(dataset,proteins,genes,conditions,colorValues,replicates,group))
}