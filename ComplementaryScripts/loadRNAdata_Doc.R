loadRNAdata_Doc <- function(dataPath,organism){
  setwd(dataPath)
  filename <- paste('STAR_',organism,'.csv',sep='')
  dataset  <-read.csv(filename,row.names = 1)
  # Add the grouping information (and enforce the order of the samples).
    conditions <- c('Ref','Doc')
    colorValues <- c("blue", "light green")
    replicates <- c(3,3)
    group <- factor(c(rep(conditions[1],replicates[1]),
                      rep(conditions[2],replicates[2])), levels = conditions)
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