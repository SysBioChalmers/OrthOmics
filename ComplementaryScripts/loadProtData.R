loadProtData <- function(dataPath,organism,emPAI){
  setwd(dataPath)
  filename <- paste(organism,'_relative_emPAI.txt',sep='')
  dataset <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  # Add the grouping information (and enforce the order of the samples).
  if (all(organism == 'yli')){
    conditions <- c('Ref','HiT','LpH')
    colorValues <- c("black", "red", "#009E73")
    replicates <- c(4,3,4)
    group <- factor(c(rep(conditions[1],replicates[1]),
                      rep(conditions[2],replicates[2]),
                      rep(conditions[3],replicates[3])), levels = conditions)
    #Get proteins and genes IDs
    IDs      <- dataset[,1:4]
    dataset  <- dataset[,(ncol(IDs)+1):(sum(replicates)+ncol(IDs))]
    dataset[is.na(dataset)] <- 0 
    genes    <- IDs[,2]
    proteins <- IDs[,4]
    proteins[proteins==''] <- genes[proteins=='']
  } else {
    conditions <- c('Ref','HiT','LpH','Osm')
    colorValues <- c("black", "red", "#009E73","blue")
    replicates <- c(3,3,3,3)
    group <- factor(c(rep(conditions[1],replicates[1]),
                      rep(conditions[2],replicates[2]),
                      rep(conditions[3],replicates[3]),
                      rep(conditions[4],replicates[4])), levels = conditions)
    #Get proteins and genes IDs
    if (all(organism =='sce')){
      IDs      <- dataset[,1:3]
      dataset  <- dataset[,(ncol(IDs)+1):(sum(replicates)+ncol(IDs))]
      dataset[is.na(dataset)] <- 0 
      genes    <- IDs[,2]
      genes[genes==''] <- IDs[genes=='',1]
      proteins <- IDs[,3]
      proteins[proteins==''] <- genes[proteins=='']
    }else{
      IDs      <- dataset[,1:3]
      dataset  <- dataset[,(ncol(IDs)+1):(sum(replicates)+ncol(IDs))]
      dataset[is.na(dataset)] <- 0 
      genes    <- IDs[,1]
      genes[genes==''] <- IDs[genes=='',1]
      proteins <- IDs[,2]
      proteins[proteins==''] <- genes[proteins=='']
    }
  }
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