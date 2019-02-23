filterLowReads <- function(data,coverage,filterType,grouping){
#Remove low reads
keep.exprs <- c()
cpm_       <- cpm(data, log = T)
for (i in 1:nrow(data)){
  #Get the 
  rowData <- cpm_[i,]
  if (all(filterType=='ref')){
    #Keep those rows which show a cpm>=1 for the std condition (on average)
    if (mean(rowData[1:grouping[1]])>=0){keep.exprs <- c(keep.exprs,i)}    
  } else{
    #Keep those rows which show an average cpm>=1 in at least (coverage) of the total samples 
    #in which it was measured
    if (length(rowData[rowData>=0])>=coverage*length(rowData)) {keep.exprs <- c(keep.exprs,i)}
   
  }
  
}
data <- data[keep.exprs,]
return(data)
}