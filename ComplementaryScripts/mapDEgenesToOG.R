mapDEgenesToOG <- function(organisms,pVal,logFC,adjustedPval,omics,repoPath){
#proteomics_Analysis
#
#Function that identifies the DE genes/proteins for a given organism in all 
#of the experimental conditions for which a DE results file is available and
#maps the genes into a OG list file (1:1:1 orthologous genes files from orthoFinder)
#
# organisms        (vector) Organism IDs (c(sce, kma or yli))
# pVal             (double) DE pValue threshold
# logFC            (double) abs(log2FC) for the DE fold-change threshold
# adjustedPval     TRUE if adjusted pValue computation should be used
# omics            (string) 'RNA' or 'proteins'
# repoPath         Main repository directory
#
# Usage: mapDEgenesToOG(organisms,pVal,logFC,adjustedPval,omics,repoPath)
#
# Last modified: Ivan Domenzain. 2019-11-27
#
  
orgColors  <- c('blue','red','yellow')
conditions <- c('HiT','LpH','Osm')
setwd(paste(repoPath,'/ComplementaryScripts',sep=''))
source('plotVennDiagram.R')
#Load OG list 
dataPath <- paste(repoPath,'/Databases',sep='')
if (all(omics=='RNA')){
  setwd(paste(repoPath,'/RNA-seq',sep=''))
  dir.create('All_organisms')
  setwd(paste(repoPath,'/RNA-seq/All_organisms',sep=''))
  resultsPath <- paste(repoPath,'/RNA-seq/All_organisms/DE_log2FC_',logFC,'_FDR_',pVal,sep='')
  dir.create(resultsPath)
} else {
  organisms  <- c('sce','kma','yli')
  setwd(paste(repoPath,'/Proteomics/Results',sep=''))
  dir.create('All_organisms')
  setwd(paste(repoPath,'/Proteomics/Results/All_organisms',sep=''))
  resultsPath <- paste(repoPath,'/Proteomics/Results/All_organisms/DE_log2FC_',logFC,'_FDR_',pVal,sep='')
  dir.create(resultsPath)
}
setwd(dataPath)
#dir.create(resultsPath)
OGlist  <- read.csv('SingleCopyOG_All.txt', header = TRUE, sep = "\t",stringsAsFactors=FALSE)

for (i in 1:length(conditions)){
  cond    <- conditions[i]
  cat(paste('Mapping DE ',omics,' for ', cond,'\n',sep=""))
  orgs    <- organisms
  DEdata  <- list()
  upReg   <- list()
  DownReg <- list()
  OGup    <- list()
  OGDown  <- list()
  if (all(cond=='Osm')){
    orgs <- organisms[1:2]
    colorValues <- orgColors[1:2]
    intLabSize <- c(rep(2.5,3))
    intLabSize[2]<-3.5
    ellipses <- 2
  }else{
    colorValues <- orgColors
    intLabSize <- c(rep(3,7))
    intLabSize[2]<-3.5
    intLabSize[4]<-3.5
    intLabSize[6]<-3.5
    intLabSize[5]<-4.5
    ellipses <- 3
  }
  for (j in 1:length(orgs)){
    print(orgs[j])
    #Load DE results file 
    if (all(omics=='RNA')){
      dataPath <- paste(repoPath,'/RNA-seq/',orgs[j],'/Results/DE_log2FC_',logFC,'_FDR_0.01',sep='')
    } else{
      dataPath <- paste(repoPath,'/Proteomics/Results/',orgs[j],'/DE_log2FC_',logFC,'_FDR_0.01',sep='')
    }
    setwd(dataPath)
    filename <- paste(orgs[j],'_',omics,'_ref_',cond,'.csv',sep='')
    #For j-th organism get the data for the genes that were DE exclusively in the i-th condition
    DEdata[[j]]  <- read.csv(filename,row.names = 1,stringsAsFactors = FALSE)
    if (adjustedPval == TRUE){
      upReg[[j]]   <- rownames(DEdata[[j]])[(DEdata[[j]]$logFC>=logFC) & (DEdata[[j]]$FDR<=pVal)]
      DownReg[[j]] <- rownames(DEdata[[j]])[(DEdata[[j]]$logFC<=-logFC) & (DEdata[[j]]$FDR<=pVal)]
    } else {
      upReg[[j]]   <- rownames(DEdata[[j]])[(DEdata[[j]]$logFC>=logFC) & (DEdata[[j]]$PValue<=pVal)]
      DownReg[[j]] <- rownames(DEdata[[j]])[(DEdata[[j]]$logFC<=-logFC) & (DEdata[[j]]$PValue<=pVal)]   
    }
    #Map the DE genes to the OG list
    k <- j+1
    #upregulated
    indexesOG  <- which(is.element(OGlist[,k],upReg[[j]]))
    indexes_up <- which(is.element(upReg[[j]],OGlist[,k]))
    OGup[[j]]  <- OGlist[indexesOG,1] 
    #Downregulated
    indexesOG     <- which(is.element(OGlist[,k],DownReg[[j]]))
    indexes_down  <- which(is.element(DownReg[[j]],OGlist[,k]))
    OGDown[[j]]   <- OGlist[indexesOG,1]  
    #How many of the DE for the j-th organism and the i-th condition were actually mapped to the OG list?
    temp <- list(upReg[[j]],upReg[[j]][indexes_up])
    #Get and save Venn diagrams for the mapping results for each organism
    setwd(resultsPath)
    png(paste(orgs[j],'_mapped_',cond,'_Up.png',sep=''),width = 600, height = 600)
    OG_mapped <- plotVennDiagram(temp,c(paste(orgs[j],'_',cond,sep=''),'OG'),c(colorValues[j],'cyan'),c(3,4,3),2)
    dev.off()
    temp <- list(DownReg[[j]],DownReg[[j]][indexes_down])
    png(paste(orgs[j],'_mapped_',cond,'_Down.png',sep=''),width = 600, height = 600)
    OG_mapped <- plotVennDiagram(temp,c(paste(orgs[j],'_',cond,sep=''),'OG'),c(colorValues[j],'cyan'),c(3,4,3),2)
    dev.off()
  }
  
  setwd(resultsPath)
  #get venn diagrams across organisms
  png(paste(omics,'_',cond,'_Up_OG.png',sep=''),width = 600, height = 600)
  conds_Up_subsets <- plotVennDiagram(OGup,orgs,colorValues,intLabSize,ellipses,TRUE)
  dev.off()
  png(paste(omics,'_',cond,'_down_OG.png',sep=''),width = 600, height = 600)
  conds_down_subsets <- plotVennDiagram(OGDown,orgs,colorValues,intLabSize,ellipses,TRUE)
  dev.off()
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
    filename <- paste('DE_',direction,'_',cond,'_OG_allOrgs.csv',sep='')
    write.table(dataCondition, filename, sep=",",row.names = FALSE,quote=FALSE)
  }
} 
cat(paste('Mapping completed, results are stored in ',omics,'/DE_log2FC_',logFC,'_FDR_',pVal,'\n',sep=""))
}

