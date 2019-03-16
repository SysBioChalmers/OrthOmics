DEpairwiseAnalysis <- function(dataset,org,conditions,coloVals,logPval,log2FC,adjusted,omics,groups,normMethod){
  nargin <- length(as.list(match.call())) -1
  if (nargin < 10){normMethod <- '-'}
  #Create new directory for storing results for the given set of threshold values 
  newDirName <- paste('DE_log2FC_',log2FC,'_FDR_',10^(-logPval),sep='')
  dir.create(newDirName)
  setwd(newDirName)

  #Group all the genes that are upregulated or downregulated in at least one of the conditions
  upReg_AllConds   <- list()
  downReg_AllConds <- list()
  DEgenes          <- c()
  data <- list()
  for (i in 2:length(conditions)) {
    print(paste('DE analysis Ref vs. ',conditions[i],sep=''))
    #Prepate DGEList object for DE analysis
    #For correct dispersion estimations it is needed to extract the just data for the reference and i-th 
    #condition
    tempConds             <- c(conditions[1],conditions[i]) # Reference first!
    colIndexes            <- groups=='Ref' | groups==conditions[i]
    dataSet               <- DGEList(counts = (dataset[,colIndexes]), genes = rownames(dataset))
    dataSet$samples$group <- groups[colIndexes]
    #Normalization as an optional step [TMM normalization is recommended for RNAseq data]
    if (all(normMethod!='-')){dataSet <- calcNormFactors(dataSet, method = normMethod)}
    dataSet               <- estimateDisp(dataSet)
    #Differential expression analysis function 
    de <- exactTest(dataSet, pair = tempConds) 
    tt <- topTags(de, n = Inf)
    # Write CSV file with the DE analysis results for all genes
    filename <- paste(org,'_',omics,'_ref_',conditions[i],'.csv',sep='')
    write.csv(tt, file = filename, row.names = F,quote = FALSE)
    # Make volcano plot highlighting significant DE hits
    if (adjusted == TRUE){
    volcanoData <- cbind(tt$table$logFC, -log10(tt$table$FDR))
    }else{volcanoData <- cbind(tt$table$logFC, -log10(tt$table$PValue))}
    colnames(volcanoData) <- c("logFC", "-log10Pval")
    rownames(volcanoData) <- tt$table$genes
    volcanoData <- as.data.frame(volcanoData)
    #Identify the DE genes from the whole list
    significantDE           <- ((abs(volcanoData$logFC)>=log2FC & volcanoData$`-log10Pval`>=logPval))
    print(paste(sum(significantDE),' DE hits',sep=''))
    #Sort DE hits by directionality of their FC
    upReg_AllConds[[i-1]]   <- rownames(volcanoData)[(volcanoData$logFC>=log2FC & volcanoData$`-log10Pval`>=logPval)]
    print(paste(length(upReg_AllConds[[i-1]]),' UpRegulated ',omics,sep=''))
    downReg_AllConds[[i-1]] <- rownames(volcanoData)[(volcanoData$logFC<=-log2FC & volcanoData$`-log10Pval`>=logPval)]
    print(paste(length(downReg_AllConds[[i-1]]),' DownRegulated ',omics,sep=''))
    #Black color for all genes
    genesColor  <- c(rep('black',nrow(dataSet)))
    #Red for downregulated
    genesColor[volcanoData$logFC<=-log2FC & volcanoData$`-log10Pval`>=logPval] = 'red'
    #Green for upregulated genes
    genesColor[volcanoData$logFC>=log2FC & volcanoData$`-log10Pval`>=logPval] = 'green'
    #Add labels to DE genes
    dataLabels <- vector(mode = "list", length = dim(volcanoData)[1])
    dataLabels[significantDE] <-tt$table$genes[significantDE]
    dataLabels[!significantDE]<-''
    DEgenes[[i-1]] <- rownames(volcanoData)[significantDE]
    yLimit <-max(volcanoData$`-log10Pval`)
    p <- ggplot(volcanoData,aes(logFC,`-log10Pval`)) 
    p <- p + geom_hline(aes(yintercept=logPval)) 
    p <- p + geom_vline(aes(xintercept=-log2FC)) + geom_vline(aes(xintercept=log2FC))
    p <- p + geom_point(colour = genesColor, alpha = 0.5)
    #p <- p + geom_text(aes(label=dataLabels),check_overlap = TRUE)
    textStr <- paste(length(downReg_AllConds[[i-1]]),' Down regulated',sep='')
    p <- p + annotate(geom = 'text',x=-0.5*log2FC, y=0.6*yLimit, label=textStr,size = 5)
    textStr <- paste(length(upReg_AllConds[[i-1]]),' Up regulated',sep='')
    p <- p + annotate(geom = 'text',x=1.5*log2FC, y=0.8*yLimit, label=textStr,size = 5)
    p <- p + labs(x='Log2FC', y = '-Log10 P-value')
    p <- p + theme_bw()
    plot(p)
    plot_name <- paste(org,'_',omics,'_ref_',conditions[i],'.png',sep='')
    ggsave(plot_name, width = 5, height=5)
    data[[i-1]] <- volcanoData
    dev.off()
  }
  #Get venn diagrams for the up and down regulated genes
  ellipses <- length(conditions)-1
  if (ellipses == 3) {
    intLabSize <- c(rep(2.5,7))
    intLabSize[2]<-3
    intLabSize[4]<-3
    intLabSize[6]<-3
    intLabSize[5]<-3.5
  }else{
    intLabSize <- c(rep(2.5,3))
    intLabSize[2]<-3.5
  }
  if (ellipses == 1) {intLabSize = 3.5}
  #Venn for upregulated genes
  png(paste(org,'_',omics,'_UpRegulated_vennAllconds.png',sep=''),width = 600, height = 600)
  allCondsUp <- plotVennDiagram(upReg_AllConds,conditions[2:length(conditions)],coloVals[2:length(conditions)],intLabSize,ellipses)
  dev.off()
  #Venn for downregulated genes
  png(paste(org,'_',omics,'_downRegulated_vennAllconds.png',sep=''),width = 600, height = 600)
  allCondsDown <- plotVennDiagram(downReg_AllConds,conditions[2:length(conditions)],coloVals[2:length(conditions)],intLabSize,ellipses)
  dev.off()
  #Put together all unique DE, up and down regulated genes for any condition
  tempDEgenes <- c()
  tempDown    <- c()
  tempUp      <- c()
  Excsv_Up    <- list()
  Excsv_down  <- list()
  condIndxes  <- c(1:(length(conditions)-1))
  #Also get a file for the DE genes in all conditions
  allCondsUp   <- allCondsUp[[ellipses+1]]
  allCondsDown <- allCondsDown[[ellipses+1]]
  allUpCols    <- c()
  allDownCols  <- c()
  for (i in 1:(length(conditions)-1)){
    #Get logFC for the core stress responses
    allUpCols   <- cbind(allUpCols,data[[i]]$logFC[is.element(rownames(data[[i]]),allCondsUp)])
    allDownCols <- cbind(allDownCols,data[[i]]$logFC[is.element(rownames(data[[i]]),allCondsDown)])
    #Get the stress responses that are exclusively for each of the conditions
    Excsv_Up[[i]]   <- upReg_AllConds[[i]]
    Excsv_down[[i]] <- downReg_AllConds[[i]]
    tempDEgenes     <- c(tempDEgenes,DEgenes[[i]])
    tempDown        <- c(tempDown,downReg_AllConds[[i]])
    tempUp          <- c(tempUp,upReg_AllConds[[i]])
    otherIndxs      <- setdiff(condIndxes,i)
    for (j in otherIndxs){
      Excsv_Up[[i]]   <- setdiff(Excsv_Up[[i]],upReg_AllConds[[j]])
      Excsv_down[[i]] <-  setdiff(Excsv_down[[i]],downReg_AllConds[[j]])
    }
    #Write a file with the genes data exclusively DE for the i-th condition
    temp <- c(Excsv_Up[[i]],Excsv_down[[i]])
    temp <- data[[i]][is.element(rownames(data[[i]]),temp),]
    filename <- paste(org,'_',omics,'_DE_exclusive_',conditions[i+1],'.csv',sep='')
    write.csv(temp, file = filename, row.names = TRUE, quote= FALSE)
  }
  #Write a file with the core stress responses (DE for all conditions)
  allUpCols <- as.data.frame(cbind(allCondsUp,allUpCols))
  filename <- paste(org,'_',omics,'_UpReg_allConds.csv',sep='')
  colnames(allUpCols) <- c('genes',conditions[2:length(conditions)])
  write.csv(allUpCols, file = filename, row.names = FALSE,quote = FALSE)

  allDownCols <- as.data.frame(cbind(allCondsDown,allDownCols))
  filename <- paste(org,'_',omics,'_DownReg_allConds.csv',sep='')
  colnames(allDownCols) <- c('genes',conditions[2:length(conditions)])
  write.csv(allDownCols, file = filename, row.names = FALSE,quote = FALSE)
  
  #Create a file with all the DE genes (for any condition) including its direction
  tempDEgenes <- unique(tempDEgenes)
  tempDown    <- unique(tempDown)
  tempUp      <- unique(tempUp)
  DEdata <- c()
  DEdata$genes <- tempDEgenes
  DEdata$direction <- c(rep('',length(tempDEgenes)))
  DEdata$direction[is.element(DEdata$genes,tempDown)] <- 'down'
  DEdata$direction[is.element(DEdata$genes,tempUp)] <- 'up'
  DEdata$direction[is.element(DEdata$genes,tempUp) & is.element(DEdata$genes,tempDown)] <- 'mixed'
  DEdata <- as.data.frame(DEdata)
  filename <- paste(org,'_',omics,'_DE_anyCondition.csv',sep='')
  write.csv(DEdata, file = filename, row.names = FALSE, quote= FALSE)
  #Create a file with the promiscuous DE genes
  promiscuous <- tempDEgenes
  promDown    <- tempDown
  promUp      <- tempUp
  for (i in 1:(length(conditions)-1)){
    promDown <- promDown[!is.element(promDown,Excsv_down[[i]])]
    promUp   <- promUp[!is.element(promUp,Excsv_Up[[i]])]
  }
  filename <- paste(org,'_',omics,'_Down_DE_promiscuous.csv',sep='')
  write.csv(promDown, file = filename, row.names = FALSE, quote= FALSE)
  filename <- paste(org,'_',omics,'_Up_DE_promiscuous.csv',sep='')
  write.csv(promUp, file = filename, row.names = FALSE, quote= FALSE)
  return(list(upReg_AllConds,downReg_AllConds,Excsv_Up,Excsv_down))
}