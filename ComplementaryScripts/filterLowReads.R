filterLowReads <- function(filtered.data,lcpm2){
#Remove low reads
keep.exprs <- c()
for (i in 1:nrow(filtered.data)){
  #Get the 
  rowData <- lcpm2[i,][!is.na(lcpm2[i,])]
  #Keep those rows which show a log2cpm >0 in at least 2/3 of the total samples 
  #in which it was measured
  if (length(rowData[rowData>0])>=(2/3)*length(rowData)){keep.exprs <- c(keep.exprs,i)}
}
filtered.data <- filtered.data[keep.exprs,]
return(filtered.data)
}