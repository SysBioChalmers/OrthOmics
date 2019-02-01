filterData <- function(dataset,grouping,metric,measured){
#Function that filters a biological dataset for different replicates in 
#different conditions. Those elements that were not measured for at least coverage
#of the triplicates for at least one of the conditions are removed. The elements
#should show low variability for all the conditions in which they were measured 
#in order to be kept
#
# Ivan Domenzain. created 2018-10-18
#
coverage <- 2/3
#First remove rows with missing IDs
NaNs    <- !is.na(rownames(dataset))
dataset <- dataset[NaNs,]
#Convert missing values to zero
dataset[is.na(dataset)] <- 0
#Save subsets for each condition
tempData <- dataset
condData <- c()
detected <- c()
for (i in 1:length(grouping)){
  condData[[i]]  <- tempData[,1:grouping[i]]
  #Identify those genes that were detected in at least coverage of the 
  #triplicates for the i-th condition
  condMat <- as.matrix(condData[[i]])
  rownames(condMat) <- c()
  detected[[i]]     <- rowSums(1*(condMat>0))
  detected[[i]]     <- which(detected[[i]]>= ((coverage)*grouping[i]))
  if (i<length(grouping)){tempData <- tempData[,(grouping[i]+1):ncol(tempData)]}
}
#Loop through all the elements (rows)
filtered <- c()
for (i in 1:nrow(dataset)){
  #print(i)
  presence  <- c(rep(FALSE,length(grouping)))
  spreading <- c(rep(FALSE,length(grouping)))
  for (j in 1:length(grouping)){
    #print(j)
    #The element should be measured in at least coverage of the replicates 
    #for being considered as present in one condition
    rowCond <- condData[[j]][i,]
    rowCond[is.na(rowCond)] <- 0
    if (sum(rowCond)>0){
      if (sum(1*(rowCond==0))<=(1/3)){presence[j] <- TRUE}
      #The element should have an RSD lower than 1 across triplicates
      #for being considered as a consistently measured value
      width <- 10
      if (all(metric == 'mean')){width <- (sd(rowCond)/mean(as.matrix(rowCond)))}
      if (all(metric == 'median')){(width <- sd(rowCond)/median(as.matrix(rowCond)))}
      if (width<=1 & width>0){spreading[j] <- TRUE}
    }
  }
  #The element should be present in at least one condition and have a low variability
  #for all the conditions in which it is present
  
  #If proteomics
  if (all(measured=='RNA')){
    conditional <- ((presence[1]==TRUE) & 1*sum(presence==TRUE) >= 0.8*length(grouping) & (all(spreading == presence)))
  }else{
    conditional <- ((presence[1]==TRUE) &  (all(spreading == presence)))
    #conditional <- (all(presence==TRUE) &  (all(spreading == presence)))
  }

  #if ((presence[1]==TRUE) & 1*sum(presence==TRUE) >= 0.8*length(grouping) & (all(spreading == presence))){
  if (conditional == TRUE){filtered <- c(filtered,i)}
}
return(list(filtered,detected))
}

