loadRNAdata <- function(dataPath,organism){
  path <- paste(dataPath,'/',organism, sep = '')
  setwd(path)
  filename <- paste(organism,'_RNA_RAW_counts.csv',sep='')
  dataset  <-read.csv(filename,row.names = 1)
  # Add the grouping information (and enforce the order of the samples).
  if (all(organism == 'yli')){
    conditions <- c('Ref','HiT','LpH')
    colorValues <- c("black", "red", "#009E73")
    replicates <- c(4,3,4)
    group <- factor(c(rep(conditions[1],replicates[1]),
                      rep(conditions[2],replicates[2]),
                      rep(conditions[3],replicates[3])), levels = conditions)
  } else {
    if (all(organism == 'cpk')){
    conditions <- c('Ref','HiT','LpH','Osm')
    colorValues <- c("black", "red", "#009E73","blue")
    replicates <- c(3,3,3,3)
    }
    if (all(organism == 'kma')){replicates <- c(3,2,3,3)}
    group <- factor(c(rep(conditions[1],replicates[1]),
                      rep(conditions[2],replicates[2]),
                      rep(conditions[3],replicates[3]),
                      rep(conditions[4],replicates[4])), levels = conditions)
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
  return(list(dataset,conditions,colorValues,replicates,group))
}