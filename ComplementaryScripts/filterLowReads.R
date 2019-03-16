filterLowReads <- function(lcpm,coverage,filterType,grouping,threshold){
nargin <- length(as.list(match.call())) -1
if (nargin < 5){threshold <- 0}
#Remove low reads
keep.exprs <- c()
for (i in 1:nrow(data)){
  #Get the 
  rowData <- lcpm[i,]
  if (all(filterType=='ref')){
    #Keep those rows which show a cpm>=1 for the std condition (on average)
    if (mean(rowData[1:grouping[1]])>=threshold){keep.exprs <- c(keep.exprs,i)}    
  } else{
    #Keep those rows which show an average cpm>=1 in at least (coverage) of the total samples 
    #in which it was measured
    if (length(rowData[rowData>=threshold])>=coverage*length(rowData)) {keep.exprs <- c(keep.exprs,i)}
  }
}
#data <- data[keep.exprs,]
return(keep.exprs)
}