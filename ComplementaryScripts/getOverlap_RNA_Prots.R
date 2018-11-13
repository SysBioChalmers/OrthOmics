repoPath  <- '/Users/ivand/Documents/GitHub/CHASSY_multiOmics_Analysis'
#Internal functions
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
setwd(scriptsPath)
source('plotVennDiagram.R')
#Provide organism code [Sce,Kma,Yli]
organism    <- 'kma'
RNAPath     <- paste(repoPath,'/RNA-seq',sep='')
protPath    <- paste(repoPath,'/Proteomics/Relative',sep='')
log2FC <- 0.5
Pvalue <- 0.05
if (all(organism == 'yli')){
  conditions <- c('Ref','HiT','LpH')
  colorValues <- c("black", "red", "#009E73")
} else {
  conditions <- c('Ref','HiT','LpH','Osm')
  colorValues <- c("black", "red", "#009E73","blue")
}

for (i in 2:length(conditions)){
  #Load DE RNA for the i-th condition
  resultsPath <- paste(RNAPath,'/',organism,'/Results',sep='')
  setwd(resultsPath)
  filename <- paste(organism,'_RNA_DE_exclusive_',conditions[i],'.csv',sep='')
  RNA <- read.delim(filename, header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
  downRNA <- RNA[which(RNA$logFC<0),]
  upRNA <- RNA[which(RNA$logFC>0),]
  #Load DE prots for the i-th condition
  resultsPath <- paste(protPath,'/Results/',organism,sep='')
  setwd(resultsPath)
  filename <- paste(organism,'_Proteins_ref_',conditions[i],'.csv',sep='')
  Prots <- read.delim(filename, header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
  downProts <- RNA[which(Prots$logFC<(-log2FC) & Prots$PValue<Pvalue),]
  upProts <- RNA[which(Prots$logFC>log2FC & Prots$PValue<Pvalue),]
  
  png(paste(organism,'_DE_RNA_prot_down_',conditions[i],'.png',sep=''),width = 600, height = 600)
  x <- plotVennDiagram(list(downRNA$X,downProts$X),c('RNA','Prots'),c(colorValues[i],'gray'),c(2.5,2.5,3.5),2)
  dev.off()
  
  png(paste(organism,'_DE_RNA_prot_up_',conditions[i],'.png',sep=''),width = 600, height = 600)
  y <- plotVennDiagram(list(upRNA$X,upProts$X),c('RNA','Prots'),c(colorValues[i],'gray'),c(2.5,2.5,3.5),2)
  dev.off()
  
}
