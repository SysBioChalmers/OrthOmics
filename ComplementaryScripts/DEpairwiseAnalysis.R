DEpairwiseAnalysis <- function(x2,org,conditions,coloVals,logPval,log2FC){
  DEgenes <- c()
  #Group all the genes that are upregulated in at least one of the conditions
  upReg_AllConds   <- c()
  #Group all the genes that are upregulated in at least one of the conditions
  downReg_AllConds <- c()
  for (i in 2:length(conditions)) {
    #Differential expression analysis
    de <- exactTest(x2, pair = c(conditions[1],conditions[i])) # Reference first!
    tt <- topTags(de, n = Inf)
    # Write CSV file
    filename <- paste(org,'_RNA_ref_',conditions[i],'.csv')
    write.csv(tt, file = filename, row.names = T)
    # Make volcano plot
    volcanoData <- cbind(tt$table$logFC, -log10(tt$table$PValue))
    colnames(volcanoData) <- c("logFC", "-log10Pval")
    rownames(volcanoData) <- tt$table$genes
    volcanoData <- as.data.frame(volcanoData)
    #Identify the DE genes from the whole list
    significantDE           <- ((abs(volcanoData$logFC)>=log2FC & volcanoData$`-log10Pval`>=logPval))
    upReg_AllConds[[i-1]]   <- rownames(volcanoData)[(volcanoData$logFC>=log2FC & volcanoData$`-log10Pval`>=logPval)]
    downReg_AllConds[[i-1]] <- rownames(volcanoData)[(volcanoData$logFC<=-log2FC & volcanoData$`-log10Pval`>=logPval)]
    
    #Black color for all genes
    genesColor  <- c(rep('black',nrow(x2)))
    #Red for downregulated
    genesColor[volcanoData$logFC<=-log2FC & volcanoData$`-log10Pval`>=logPval] = 'red'
    #Green for upregulated genes
    genesColor[volcanoData$logFC>=log2FC & volcanoData$`-log10Pval`>=logPval] = 'green'
    #Add labels to DE genes
    dataLabels <- vector(mode = "list", length = dim(volcanoData)[1])
    dataLabels[significantDE] <-tt$table$genes[significantDE]
    dataLabels[!significantDE]<-''
    DEgenes[[i-1]] <- rownames(volcanoData)[significantDE]
    
    p <- ggplot(volcanoData,aes(logFC,`-log10Pval`)) 
    p <- p + geom_hline(aes(yintercept=logPval)) 
    p <- p + geom_vline(aes(xintercept=-log2FC)) + geom_vline(aes(xintercept=log2FC))
    p <- p + geom_point(colour = genesColor, alpha = 0.5)
    #p <- p + geom_text(aes(label=dataLabels),check_overlap = TRUE)
    textStr <- paste(length(downReg_AllConds[[i-1]]),' Down regulated',sep='')
    p <- p + annotate(geom = 'text',x=-1, y=2.5, label=textStr,size = 5)
    textStr <- paste(length(upReg_AllConds[[i-1]]),' Up regulated',sep='')
    p <- p + annotate(geom = 'text',x=1, y=2.8, label=textStr,size = 5)
    p <- p + labs(x='Log2FC', y = '-Log10 P-value')
    p <- p + theme_bw()
    plot(p)
    plot_name <- paste(org,'_RNA_ref_',conditions[i],'.png',sep='')
    ggsave(plot_name, width = 5, height=5)
    dev.off()
  }
  #Get venn diagrams for the up and down regulated genes
  if (all(org != 'yli')) {
    intLabSize <- c(rep(2.5,7))
    intLabSize[2]<-3
    intLabSize[4]<-3
    intLabSize[6]<-3
    intLabSize[5]<-3.5
    ellipses <- 3
  }else{
    intLabSize <- c(rep(2.5,3))
    intLabSize[2]<-3.5
    ellipses <- 2
  }
  #Venn for upregulated genes
  png(paste(org,'_RNA_UpRegulated_vennAllconds.png',sep=''),width = 600, height = 600)
  allCondsUp <- plotVennDiagram(upReg_AllConds,conditions[2:length(conditions)],coloVals[2:length(conditions)],intLabSize,ellipses)
  dev.off()
  #Venn for downregulated genes
  png(paste(org,'_RNA_downRegulated_vennAllconds.png',sep=''),width = 600, height = 600)
  allCondsDown <- plotVennDiagram(downReg_AllConds,conditions[2:length(conditions)],coloVals[2:length(conditions)],intLabSize,ellipses)
  dev.off()
  
  #Put together all unique DE, up and down regulated genes for any condition
  tempDEgenes <- c()
  tempDown <- c()
  tempUp <- c()
  for (i in 1:(length(conditions)-1)){
    tempDEgenes <- c(tempDEgenes,DEgenes[[i]])
    tempDown <- c(tempDown,downReg_AllConds[[i]])
    tempUp <- c(tempUp,upReg_AllConds[[i]])
  }
  tempDEgenes <- unique(tempDEgenes)
  tempDown <- unique(tempDown)
  tempUp <- unique(tempUp)
  #Create a file with all the DE genes (for any condition) including its direction
  DEdata <- c()
  DEdata$genes <- tempDEgenes
  DEdata$direction <- c(rep('',length(tempDEgenes)))
  DEdata$direction[is.element(DEdata$genes,tempDown)] <- 'down'
  DEdata$direction[is.element(DEdata$genes,tempUp)] <- 'up'
  DEdata$direction[is.element(DEdata$genes,tempUp) & is.element(DEdata$genes,tempDown)] <- 'mixed'
  DEdata <- as.data.frame(DEdata)
  filename <- paste(org,'_RNA_DE_anyCondition.csv')
  write.csv(DEdata, file = filename, row.names = FALSE)
  return(list(upReg_AllConds,downReg_AllConds))
}