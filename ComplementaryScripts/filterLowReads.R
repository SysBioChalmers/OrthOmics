filterLowReads <- function(filtered.data,cpm_,Pmethod){
#Remove low reads
keep.exprs <- c()
#if (Pmethod==1){factor <- 1}
#else{
#  factor <- 2/3
#  }
factor <- 1
for (i in 1:nrow(filtered.data)){
  #Get the 
  rowData <- cpm_[i,]#[!is.na(cpm_[i,])]
  #Keep those rows which show a cpm>=1 in at least (factor) of the total samples 
  #in which it was measured
  if (length(rowData[rowData>=1])>=factor*length(rowData)){keep.exprs <- c(keep.exprs,i)}
}
filtered.data <- filtered.data[keep.exprs,]
return(filtered.data)
}