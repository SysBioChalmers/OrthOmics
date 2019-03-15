getPCAplot <- function(data,conditions,group,replicates,colVals,org,plot_name,omicsName){
  data[is.na.data.frame(data)] <- 0
  t.data<-t(as.matrix(data))
  prots.pca <- prcomp(t.data, center = TRUE,scale. = TRUE) 
  summary(prots.pca)  
  #Plot PC1 and PC2 
  labelStr <- c()
  for (i in 1:length(conditions)) {
    labelStr <- c(labelStr, seq(1,replicates[i],1))
    if (all(colVals[[i]]=='grey')){colVals[[i]]='black'}
  }
  print(colVals)
  p<-ggbiplot(prots.pca, choices = c(1,2),var.axes = FALSE,groups = group,
              varname.size = 6,labels=labelStr,labels.size = 6,ellipse=TRUE)
  p <- p + theme_bw()
  #Format plot
  p <- p + scale_colour_manual(values=colVals)
  titleStr <- paste(org,'/ ',length(data[,2]), omicsName)
  p <- p + labs(title = titleStr)
  p
  #Save plot
  ggsave(plot_name, width = 5, height=5)
  return(prots.pca)
}