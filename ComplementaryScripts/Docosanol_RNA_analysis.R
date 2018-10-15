setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/ComplementaryScripts')
source('convert_strainIDs.r')
source('plotVennDiagram.r')

organisms <- c('sce','kma','yli')
#Dataset strains
strains   <- c('s288c','DMKU','CLIB122')
#Load OrthoList
setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Databases')
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
FC     <- 0.5

for (i in 1:length(organisms)) {
  #Get DE set for the organism
  org <- organisms[i]
  setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Biopetrolia/Docosanol-RNAseq')
  filename <- paste('RNAseq_',org,'_doc.txt',sep='')
  dataset  <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE)
  print(org)
  print(paste(length(dataset$genes),' genes'))
  
  #Get differentially expressed genes
  dataset <- dataset[(dataset$logFC>=FC | dataset$logFC<=-FC),] 
  dataset <- dataset[dataset$FDR<=pValue,]
  print(paste(length(dataset$genes),' DE genes'))
  #Write file with the significantly DE genes for the organism
  setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Biopetrolia/Results')
  filename <- paste('RNAseq_DE_',org,'_docosanol.txt',sep='')
  write.table(dataset, filename, sep="\t",row.names = FALSE)
  
  #convert IDs for compatibility with the OrthoList 
  if ((all(tolower(org) == 'yli')) | all(tolower(org) == 'kma')) {
    dataset     <- convert_strainIDs(org,strains[i],dataset)
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
  setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Biopetrolia/Results')
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
direction    <- 'UpReg'
colorValues  <- colorValues <- c("blue","red", "yellow")
intLabSize   <- c(rep(2,7))
intLabSize[5]<- 3
intLabSize[2]<- 2.5
intLabSize[4]<- 2.5
intLabSize[6]<- 2.5
ellipses <- length(organisms)
indexes <- plotVennDiagram(OGDown,organisms,colorValues,intLabSize,ellipses)
#Write files for the different overlaps
for (i in 1:length(organisms)){
  if (i<3){
    j <- i+1
    k <- i
    }
  else{
    j <-1
    k <- 4
    }
  matches     <- match(indexes[[k]],orthoList[,1])
  converted   <- matches[!is.na(matches)]
  data        <- orthoList[converted,]
  filename    <- paste('DE_',direction,'_OG_',organisms[i],'_',organisms[j],'_docosanol.txt',sep='')
  write.table(data, filename, sep="\t",row.names = FALSE)
}
