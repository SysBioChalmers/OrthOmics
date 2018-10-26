#Write your repo directory here!!
repoPath <- '/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis'
setwd(paste(repoPath,'/ComplementaryScripts',sep=''))
source('convert_strainIDs.r')
source('plotVennDiagram.r')

organisms <- c('sce','kma','yli')
#Dataset strains
strains   <- c('s288c','DMKU','CLIB122')
#Load OrthoList
setwd(paste(repoPath,'/Databases',sep=''))
ListFile <- 'SingleCopyOG_All.txt'
orthoList <- read.delim(ListFile, header = TRUE, sep = "\t",stringsAsFactors=FALSE)
colnames(orthoList) <- c('OG_IDs','sce','kma','yli')
#Standardize format
orthoList$kma <- gsub('KMXK','KMXK_',orthoList$kma)

#Create necessary variables
newIDs <- c()
UpRegulated   <- c()
DownRegulated <- c()
OGUp   <- c()
OGDown <- c()
#DE thresholds
pValue <- 0.05
FC     <- 0.75

for (i in 1:length(organisms)) {
  #Get DE set for the organism
  org <- organisms[i]
  setwd(paste(repoPath,'/Final_products_Tolerance/Docosanol/Results/',org,sep=''))
  filename <- paste(org,'_RNA_ref_Doc.csv',sep='')
  dataset  <- read.delim(filename, header = TRUE, sep = ",",stringsAsFactors=FALSE)
  print(org)
  print(paste(length(dataset$genes),' genes'))
  
  #Get differentially expressed genes
  dataset <- dataset[(dataset$logFC>=FC | dataset$logFC<=-FC),] 
  dataset <- dataset[dataset$PValue<=pValue,]
  print(paste(length(dataset$genes),' DE genes'))
  #Write file with the significantly DE genes for the organism
  setwd(paste(repoPath,'/Final_products_Tolerance/Docosanol/Results',sep=''))
  filename <- paste('RNAseq_DE_',org,'_docosanol.txt',sep='')
  write.table(dataset, filename, sep="\t",row.names = FALSE)
  
  #convert IDs for compatibility with the OrthoList 
  if ((all(tolower(org) == 'yli')) | all(tolower(org) == 'kma')) {
    dataset     <- convert_strainIDs(org,strains[i],dataset,repoPath)
    newIDs[[i]] <- dataset$genes
  }
  
  #Map the dataset into the orthoList
  matches     <- match(orthoList[,i+1],dataset$genes)
  converted   <- matches[!is.na(matches)]
  matchedData <- dataset[converted,]
  print(paste(length(converted), ' in OrthoList'))
  #Substitute IDs by OG 
  matches     <- match(dataset$genes,orthoList[,i+1])
  converted   <- matches[!is.na(matches)]
  temp        <- cbind(orthoList$OG_IDs[converted],matchedData)
  #Set path for results storage
  setwd(paste(repoPath,'/Final_products_Tolerance/Docosanol/Results/Allorgs',sep=''))
  filename <- paste('RNAseq_DE_',org,'_orthologs_docosanol.txt',sep='')
  write.table(temp, filename, sep="\t",row.names = FALSE)
  #Substitute gene names for OG id's
  matchedData$genes <- orthoList$OG_IDs[converted]

  #Get up and down regulated genes from the different subsets
  UpRegulated[[i]]   <- dataset$genes[dataset$logFC>=FC] 
  DownRegulated[[i]] <- dataset$genes[dataset$logFC<=FC] 
  OGUp[[i]]          <- matchedData$genes[matchedData$logFC>=FC] 
  OGDown[[i]]        <- matchedData$genes[matchedData$logFC<=FC] 
  print(paste(length(OGUp[[i]]), '  UpRegulated'))
  print(paste(length(OGDown[[i]]), '  DownRegulated'))
  print(c())
  
}

#Get Venn diagram for common stress responses
organisms    <- c('Sce','Kma','Yli')
colorValues  <- colorValues <- c("blue","red", "yellow")
intLabSize   <- c(rep(2.5,7))
intLabSize[5]<- 3.5
intLabSize[2]<- 3
intLabSize[4]<- 3
intLabSize[6]<- 3
ellipses <- length(organisms)

setwd(paste(repoPath,'/Final_products_Tolerance/Docosanol/Results/Allorgs',sep=''))
for (direction in c('DownReg','UpReg')){
  print(direction)

  if (all(direction == 'DownReg')) {
    subsetGenes <- OGDown
    figure <- 'Docosanol_Down_venn_plot.png'
  }else{
    subsetGenes <- OGUp
    figure <- 'Docosanol_Up_venn_plot.png'
    }

  png(figure, width = 500, height = 500)
  indexes <- plotVennDiagram(subsetGenes,organisms,colorValues,intLabSize,ellipses)
  dev.off()
  #Write files for the different overlaps
  for (i in 1:length(organisms)){
    k <- i
    j <- i+1
    if (i==3){j <-1}
    columns <- c(1,i+1,j+1)
    matches     <- match(indexes[[k]],orthoList[,1])
    converted   <- matches[!is.na(matches)]
    data        <- orthoList[converted,columns]
    filename    <- paste('DE_',direction,'_OG_',organisms[i],'_',organisms[j],'_docosanol.txt',sep='')
    write.table(data, filename, sep="\t",row.names = FALSE)
  }
  matches     <- match(indexes[[4]],orthoList[,1])
  converted   <- matches[!is.na(matches)]
  data        <- orthoList[converted,]
  filename    <- paste('DE_',direction,'_OG_AllOrgs_docosanol.txt',sep='')
  write.table(data, filename, sep="\t",row.names = FALSE)
}
