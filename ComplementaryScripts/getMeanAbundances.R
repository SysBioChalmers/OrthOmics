getMeanAbundances <- function(dataset,group,conditions,filename){
meanValues <- c()
for (i in 1:length(conditions)){
  columns    <- which(group==conditions[i])
  temp       <- dataset[,columns]
  meanValues <- cbind(meanValues,rowMeans(temp))
}
rownames(meanValues) <- rownames(dataset)
write.csv(meanValues, file = filename, row.names = T,quote = FALSE)
}