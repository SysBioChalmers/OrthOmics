filterLowReads <- function(data,coverage,filterType,grouping){
#Remove low reads
keep.exprs <- c()
cpm_       <- cpm(data, log = T)
for (i in 1:nrow(data)){
  #Get the 
  rowData <- cpm_[i,]#[!is.na(cpm_[i,])]
  if (all(filterType=='all')){
      #Keep those rows which show a cpm>=1 in at least (coverage) of the total samples 
      #in which it was measured
      if (length(rowData[rowData>=1])>=coverage*length(rowData)) {keep.exprs <- c(keep.exprs,i)}
  } else{
    #Keep those rows which show a cpm>=1 for at least 2/3 of the triplicates
    #in the standard condition
    if (length(rowData[1:grouping[1]]>=1)>=(2/3)*grouping[1]) {keep.exprs <- c(keep.exprs,i)}    
  }
  
}
data <- data[keep.exprs,]
return(data)
}