library("rlist")
repoPath   <- '/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis'
organisms  <- c('sce','kma','yli')
orgColors  <- c('blue','red','yellow')
setwd(paste(repoPath,'/ComplementaryScripts',sep=''))
source('plotVennDiagram.R')

#Load OG list 
dataPath <- paste(repoPath,'/Databases',sep='')
resultsPath <- paste(repoPath,'/RNA-seq/All_organisms',sep='')
setwd(dataPath)
OGlist  <- read.csv('SingleCopyOG_All.txt', header = TRUE, sep = "\t",stringsAsFactors=FALSE)
#Loop through all the organisms 
coreGenes <- list()

for (dir in c('Down','Up')){
for (i in 1:length(organisms)){
  print(organisms[i])
  dataPath <- paste(repoPath,'/RNA-seq/',organisms[i],'/Results',sep='')
  setwd(dataPath)
  filename  <- paste(organisms[i],'_RNA_',dir,'Reg_allConds.csv',sep='')
  data      <- read.csv(filename,row.names = 1,stringsAsFactors = FALSE)
  coreGenes [[i]] <- rownames(data)
  j <- i+1
  indexesOG <- which(is.element(OGlist[,j],coreGenes[[i]]))
  
  png(paste(organisms[i],'_RNAseq_',dir,'_AllConds_OG.png',sep=''),width = 600, height = 600)
  plotVennDiagram(list(coreGenes[[i]],OGlist[indexesOG,j]),c(organisms[i],'OG'),c(orgColors[i],'gray'),c(2.5,3.5,2.5),2)
  dev.off()
  coreGenes [[i]] <- OGlist[indexesOG,1]
  
}
  setwd(resultsPath)
  
  intLabSize <- c(rep(2.5,7))
  intLabSize[2]<-3
  intLabSize[4]<-3
  intLabSize[6]<-3
  intLabSize[5]<-3.5
  
  png(paste('AllOrgs_RNAseq_',dir,'_AllConds_OG.png',sep=''),width = 600, height = 600)
  plotVennDiagram(coreGenes,organisms,orgColors,intLabSize,3)
  dev.off()
}
