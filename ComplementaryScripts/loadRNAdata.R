loadRNAdata <- function(dataPath,organism,anox){
  nargin <- length(as.list(match.call())) -1
  if (nargin < 3){anox <- FALSE}
  path <- paste(dataPath,'/',organism, sep = '')
  setwd(path)
  filename <- paste(organism,'_RNA_RAW_counts.csv',sep='')
  #dataset  <-read.csv(filename,row.names = 1)
  dataset <- read.delim(filename, header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
  # Add the grouping information (and enforce the order of the samples).
  if (all(organism == 'yli')){
    conditions <- c('Ref','HiT','LpH')
    colorValues <- c("grey", "red", "#009E73")
    replicates <- c(4,3,4)
    excl <- c()
  } else {
    if (all(organism == 'sce')){
      conditions <- c('Ref','HiT','LpH','Osm','Ana')
      colorValues <- c("grey", "red", "#009E73","blue",'yellow')
      replicates <- c(3,3,3,3,3)
      if (anox){
        excl <- which(conditions!='Ana' & conditions!='Ref')
      }else{
        excl <- which(conditions=='Ana')
      }
    }
    if (all(organism == 'kma')){
      replicates  <- c(3,2,3,3)
      conditions  <- c('Ref','HiT','LpH','Osm')
      colorValues <- c("grey", "red", "#009E73","blue")
      excl <- c()
      }
  }
  
  #Generate colNames and grouping
  cols <- c()
  groupVector <- c()
  for (i in 1:length(conditions)) {
    for (j in 1:replicates[i]) {
      cols        <- c(cols,paste(conditions[i],'_',j,sep=''))
    }
    groupVector <- c(groupVector,c(rep(conditions[i],replicates[i])))
  }
  #Get counts dataset and genes
  genes   <- dataset[,1]
  dataset <- dataset[,2:ncol(dataset)]
  colnames(dataset)       <- cols
  dataset[is.na(dataset)] <- 0
  #Exclude the undesired condition columns
  for (exclCol in excl){
    #Get columns to remove in reamining dataset
    toRemove <- which(groupVector==conditions[exclCol])
    #Update dataset
    dataset     <- dataset[,-toRemove]
    groupVector <- groupVector[-toRemove]
  }
  #Update conditions, colors and replicates
  conditions  <- conditions[-excl]
  colorValues <- colorValues[-excl]
  replicates  <- replicates[-excl]
  
  group <- factor(groupVector,levels = conditions)
  remove(cols)
  return(list(dataset,genes,conditions,colorValues,replicates,group))
}