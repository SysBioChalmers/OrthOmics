repoPath  <- '/Users/ivand/Documents/GitHub/CHASSY_multiOmics_Analysis'
#Internal functions
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
setwd(scriptsPath)
source('plotVennDiagram.R')
#Provide organism code [Sce,Kma,Yli]
organism    <- 'kma'
if (all(organism=='sce')){
  org = 'cpk'
} else{org = organism}
RNAPath     <- paste(repoPath,'/RNA-seq',sep='')
protPath    <- paste(repoPath,'/Proteomics',sep='')
log2FC <- 1
Pvalue <- 0.01
adjProteomics <- TRUE
if (all(organism == 'yli')){
  conditions <- c('Ref','HiT','LpH')
  colorValues <- c("grey", "red", "#009E73")
} else {
  conditions <- c('Ref','HiT','LpH','Osm')
  colorValues <- c("grey", "red", "#009E73","blue")
}

for (i in 2:length(conditions)){
  #Load DE RNA for the i-th condition
  resultsPath <- paste(RNAPath,'/',org,'/Results/',sep='')
  resultsPath <- paste(resultsPath,'DE_log2FC_',log2FC,'_FDR_',Pvalue,sep='')
  setwd(resultsPath)
  filename <- paste(org,'_RNA_ref_',conditions[i],'.csv',sep='')
  RNA      <- read.delim(filename, header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
  downRNA  <- RNA$genes[which(RNA$logFC<(-log2FC) & RNA$FDR<Pvalue)]
  upRNA    <- RNA$genes[which(RNA$logFC>(-log2FC) & RNA$FDR<Pvalue)]
  #Load DE prots for the i-th condition
  resultsPath <- paste(protPath,'/Results/',organism,sep='')
  resultsPath <- paste(resultsPath,'/DE_log2FC_',log2FC,'_FDR_',Pvalue,sep='')
  setwd(resultsPath)
  filename <- paste(organism,'_Proteins_ref_',conditions[i],'.csv',sep='')
  Prots    <- read.delim(filename, header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
  if (all(organism == 'kma')){Prots$genes <- gsub("KMXK", "KMXK_", Prots$genes)}
  if (adjProteomics){
    downProts <- Prots$genes[which(Prots$logFC<(-log2FC) & Prots$FDR<Pvalue)]
    upProts   <- Prots$genes[which(Prots$logFC>log2FC & Prots$FDR<Pvalue)]
  } else{
    downProts <- Prots$genes[which(Prots$logFC<(-log2FC) & Prots$PValue<Pvalue)]
    upProts   <- Prots$genes[which(Prots$logFC>log2FC & Prots$PValue<Pvalue)]
  }

  png(paste(organism,'_DE_RNA_prot_down_',conditions[i],'.png',sep=''),width = 600, height = 600)
  x <- plotVennDiagram(list(downRNA,downProts),c('RNA','Prots'),c(colorValues[i],'gray'),c(3,4,3),2)
  dev.off()
  
  png(paste(organism,'_DE_RNA_prot_up_',conditions[i],'.png',sep=''),width = 600, height = 600)
  y <- plotVennDiagram(list(upRNA,upProts),c('RNA','Prots'),c(colorValues[i],'gray'),c(3,4,3),2)
  dev.off()
  
  filename <- paste(organism,'_DEgenes_prot_RNAoverlap_',conditions[i],'.csv',sep='')
  table <- as.data.frame(c(x[[3]],y[[3]]))
  colnames(table) <- paste('Ref_',conditions[i],sep='')
  write.csv(table, file = filename, row.names = FALSE)
}
