install.packages("rlist")
library("rlist")
repoPath   <- '/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis'
organisms  <- c('sce','kma','yli')
orgColors  <- c('blue','red','yellow')

conditions <- c('HiT','LpH','Osm')
setwd(paste(repoPath,'/ComplementaryScripts',sep=''))
source('plotVennDiagram.R')
#Load OG list 
dataPath <- paste(repoPath,'/Databases',sep='')
resultsPath <- paste(repoPath,'/RNA-seq/All_organisms',sep='')
setwd(dataPath)
OGlist  <- read.csv('SingleCopyOG_All.txt', header = TRUE, sep = "\t",stringsAsFactors=FALSE)
#DE thresholds
logPVal <- abs(log10(0.01))
logFC <- 0.75
for (i in 1:length(conditions)){
  cond    <- conditions[i]
  print(cond)
  orgs    <- organisms
  DEdata  <- list()
  upReg   <- list()
  DownReg <- list()
  OGup    <- list()
  OGDown  <- list()
  if (i==3){
    orgs <- organisms[1:2]
    colorValues <- orgColors[1:2]
    intLabSize <- c(rep(2.5,3))
    intLabSize[2]<-3.5
    ellipses <- 2
  }else{
    colorValues <- orgColors
    intLabSize <- c(rep(2.5,7))
    intLabSize[2]<-3
    intLabSize[4]<-3
    intLabSize[6]<-3
    intLabSize[5]<-3.5
    ellipses <- 3
  }
  for (j in 1:length(orgs)){
    print(orgs[j])
    dataPath <- paste(repoPath,'/RNA-seq/',orgs[j],'/Results',sep='')
    setwd(dataPath)
    filename     <- paste(orgs[j],'_RNA_DE_exclusive_',cond,'.csv',sep='')
    #For j-th organism get the data for the genes that were DE exclusively in the i-th condition
    DEdata[[j]]  <- read.csv(filename,row.names = 1,stringsAsFactors = FALSE)
    upReg[[j]]   <- rownames(DEdata[[j]])[(DEdata[[j]]$logFC>=logFC) & (DEdata[[j]]$X.log10Pval>=logPVal)]
    DownReg[[j]] <- rownames(DEdata[[j]])[(DEdata[[j]]$logFC<=-logFC) & (DEdata[[j]]$X.log10Pval>=logPVal)]
    #Map the DE genes to the OG list
    k <- j+1
    #upregulated
    #OGup[[j]] <- upReg[[j]]
    indexesOG <- which(is.element(OGlist[,k],upReg[[j]]))
    indexes_up   <- which(is.element(upReg[[j]],OGlist[,k]))
    #OGup[[j]][indexes] <- OGlist[indexesOG,1] 
    OGup[[j]] <- OGlist[indexesOG,1] 
    #Downregulated
    #OGDown[[j]] <- DownReg[[j]]
    indexesOG <- which(is.element(OGlist[,k],DownReg[[j]]))
    indexes_down   <- which(is.element(DownReg[[j]],OGlist[,k]))
    #OGDown[[j]][indexes] <- OGlist[indexesOG,1]  
    OGDown[[j]] <- OGlist[indexesOG,1]  
    #How many of the DE for the j-th organism and the i-th condition were actually mapped to the OG list?
    temp <- list(upReg[[j]],upReg[[j]][indexes_up])
    setwd(resultsPath)
    png(paste(orgs[j],'_mapped_',cond,'_Exclusive_Up.png',sep=''),width = 600, height = 600)
    OG_mapped <- plotVennDiagram(temp,c(paste(orgs[j],'_',cond,sep=''),'OG'),c(colorValues[j],'cyan'),c(2.5,3.5,2.5),2)
    dev.off()
    temp <- list(DownReg[[j]],DownReg[[j]][indexes_down])
    png(paste(orgs[j],'_mapped_',cond,'_Exclusive_Down.png',sep=''),width = 600, height = 600)
    OG_mapped <- plotVennDiagram(temp,c(paste(orgs[j],'_',cond,sep=''),'OG'),c(colorValues[j],'cyan'),c(2.5,3.5,2.5),2)
    dev.off()
  }
  
  setwd(resultsPath)
  #get venn diagrams across organisms
  png(paste('RNAseq_',cond,'_Exclusive_Up_OG.png',sep=''),width = 600, height = 600)
  conds_Up_subsets <- plotVennDiagram(OGup,orgs,colorValues,intLabSize,ellipses)
  dev.off()
  png(paste('RNAseq_',cond,'_Exclusive_down_OG.png',sep=''),width = 600, height = 600)
  conds_down_subsets <- plotVennDiagram(OGDown,orgs,colorValues,intLabSize,ellipses)
  dev.off
  #Write files for the different overlaps
  for (direction in c('up','down')){
    if (all(direction=='up')){
      subsets <- OGup
    }else{subsets  <- OGDown}
    #Merge all the DE OG for all the organisms on the 
    condOG_allOrgs <- list.cases(subsets)
    genesOrg       <- list()
    #Map each OG to their correspondant organism genes
    for (k in 1:length(orgs)){
      genesOrg[[k]]  <- c(rep('',length(condOG_allOrgs)))
      OGcond_indexes <- match(subsets[[k]],condOG_allOrgs)
      genesOrg[[k]][OGcond_indexes] <- OGlist[match(condOG_allOrgs[OGcond_indexes],OGlist[,1]),k+1]
    }
    dataCondition <- as.data.frame(genesOrg)
    dataCondition <- cbind(condOG_allOrgs,dataCondition)
    colnames(dataCondition) <- c('OG',orgs)
    filename <- paste('DE_',direction,'_exclusive_',cond,'_OG_allOrgs.csv',sep='')
    write.table(dataCondition, filename, sep=",",row.names = FALSE,quote=FALSE)
  }
}  


