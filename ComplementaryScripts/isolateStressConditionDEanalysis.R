#Get DE genes specific to a given condition
organisms    <- c('sce','yli','kma')
repoPath     <- '/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis'
setwd(paste(repoPath,'/ComplementaryScripts',sep=''))
source('plotVennDiagram.r')
logFC    <- 0.75
pVal     <- 0.01

condUp   <- list()
condDown <- list()
Exc_up   <- list()
Exc_Down <- list()
subsets  <- list()

for (i in 1:length(organisms)){
  #Load other stresses data
  OrgDataPath  <- paste(repoPath,'/RNA-seq/',organisms[i],'/Results',sep='')
  setwd(OrgDataPath)
  filename <- paste(organisms[i],'_RNA_DE_anyCondition.csv')
  data     <- read.csv(filename,row.names = 1,stringsAsFactors = FALSE)
  Allup    <- rownames(data)[data$direction=='up' |data$direction=='mixed']
  Alldown  <- rownames(data)[data$direction=='down' |data$direction=='mixed']
  #Load Docosanol DE data
  DataPath <- paste(repoPath,'/Final_products_Tolerance/Docosanol/Results/',organisms[i],sep='')
  setwd(DataPath)
  filename <- paste(organisms[i],'_RNA_ref_Doc.csv',sep='')
  DEdata   <- read.csv(filename,row.names = 1,stringsAsFactors = FALSE)
  condUp[[i]]   <- DEdata$genes[DEdata$logFC>=logFC & DEdata$PValue <= pVal]
  condDown[[i]] <- DEdata$genes[DEdata$logFC<=-logFC & DEdata$PValue <= pVal]
  Exc_up[[i]]   <- condUp[[i]][!is.element(condUp[[i]],Allup)] 
  Exc_Down[[i]] <- condDown[[i]][!is.element(condDown[[i]],Alldown)] 
  directions <- c('up','down')
  for (j in 1:2){
    direction <- directions[j]
    if (all(direction=='up')){elementsList = list(condUp[[i]],Allup)}
    else{elementsList = list(condDown[[i]],Alldown)}
    setwd(DataPath)
    figure  <- paste(organisms[i],'_Doc_vs_OtherStresses_',direction,'.png',sep='')
    png(figure, width = 500, height = 500)
    subsets <- plotVennDiagram(elementsList,c('Doc','Other'),c('light green','red'),c(2.5,3,2.5),2)
    dev.off()
  }
  #Write file with the exclusive stress responses to Docosanol
  data <- DEdata[(is.element(DEdata$genes,condUp[[i]]) |is.element(DEdata$genes,condDown[[i]])),]
  filename <- paste(organisms[i],'_Docosanol_specific_DEgenes.csv',sep='')
  write.csv(data, file = filename, row.names = TRUE,quote = FALSE)
}