createIntegratedTable <- function(organism,FDR,log2FC,adjustedP,justDEgenes,repoPath){
#
#Function that integrates the results from RNAseq DE with absolute proteomics
#measurements and genes annotation.
#
# Ivan Domenzain. Last modified 2019-11-26
#
  
#==================================== DEFINE VARIABLES =======================================
#Provide organism code [sce,kma,yli]
if (all(organism == 'yli')){
  conditions  <- c('HiT','LpH')
  colorValues <- c("red", "#009E73")
} else {
  conditions  <- c('HiT','LpH','Osm')
  colorValues <- c("red", "#009E73","blue")
}
#Define DE thresholds
logPval     <- abs(log10(pVal))
#=============================== Relevant directories ======================================
#Relevant paths (The user should provide the path in which the repository is stored)
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
resultsPath <- paste(repoPath,'/Integration',sep='')
uniprotPath <- paste(repoPath,'/Databases/Uniprot',sep='')
RNApath     <- paste(repoPath,'/RNA-seq',sep='')
Protpath    <- paste(repoPath,'/Proteomics/Results/',organism,sep='')
Orthpath    <- paste(repoPath,'/Orthologs',sep='')
GOTermsPath <- paste(repoPath,'/Databases/GO_terms',sep='')
#=============================== Load data =================================================
#============= Load absolute proteomics measurements (all conditions)
setwd(Protpath)
fileName <- paste(organism,'_abs_NSAF_filtered.txt',sep="")
prot_Abs <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE)
proteins <- rownames(prot_Abs)
#============= Load gene groups (orthology) list
setwd(Orthpath)
fileName   <- paste(organism,'_orthology_groups.txt',sep='')
geneGroups <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
if (all(organism=='kma')){
  proteins        <- gsub('KMXK','KMXK_',proteins)
  uniprotFile     <- paste('uniprot_',organism,'_DMKU.txt',sep='')
  #mapping between DMKU and CBS6556 gene IDs
  fileName        <- 'kma_strains_orthologs.txt'
  strainOrthologs <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
} else{
  if (all(organism=='yli')){
    uniprotFile <- paste('uniprot_',organism,'_W29.txt',sep='')
  } else {
      if (all(organism=='sce')){
        uniprotFile <- paste('uniprot_sce.txt',sep='')
      } 
  }
}
#============= Load GO terms file
setwd(GOTermsPath)
fileName <- paste(organism,'_GOterms.txt',sep="")
GO_terms <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
#============= Load uniprot file
setwd(uniprotPath)
uniprotData  <- read.delim(uniprotFile, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")

#Loop through all conditions 
cat(paste("Creating intergrated Omics table for: ",organism,'\n',sep=""))
tableData <- c()
for (i in 1:length(conditions)){
  condition <- conditions[i]
  cat(paste("Adding DE genes for ",condition,'\n',sep=""))
  cat("Adding semi-absolute proteomics data\n")
  cat("Adding gene orthologous grouping info\n")
  cat("Adding gene uniprot info\n")
  cat("Adding GO terms info\n")
  #Load RNAseq DE results for the i-th condition
  setwd(paste(RNApath,'/',organism,'/Results/DE_log2FC_',log2FC,'_FDR_',FDR,sep=''))
  fileName    <- paste(organism,'_RNA_ref_',conditions[i],'.csv',sep='')
  RNA_DE_data <- read.delim(fileName, header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
  genes       <- RNA_DE_data$genes
  for (j in 1:length(genes)){
    row          <- c()
    gene         <- genes[j]
    originalGene <- gene
    #Identify if j-th was significantly DE
    significance <- FALSE
    if (RNA_DE_data$FDR[j] < FDR & abs(RNA_DE_data$logFC[j])>=log2FC){
        significance <- TRUE
    }
    direction <- sign(RNA_DE_data$logFC[j])
    #Search protein abundances for queried gene
    ProtIndex <- which(proteins == gene)
    if (length(ProtIndex)>0){
      P_ref  <- prot_Abs[ProtIndex,1]
      P_cond <- prot_Abs[ProtIndex,i+1]
    } else {
      P_ref  <- NaN
      P_cond <- NaN
    }
    #Search gene orthoGroup
    geneIndex <- which(geneGroups[,1] == gene)
    orthoGene <- c()
    if (length(geneIndex)>0){
      orthoGroup <- geneGroups[geneIndex,2]
      #Map CBS6556 IDs to DMKU 
      if (all(organism == 'kma')){
          geneIndex    <- which(strainOrthologs[,2] == gene)[1]
          originalGene <- gene
          gene         <- strainOrthologs[geneIndex,1]
          orthoGene    <- originalGene
      } else {orthoGene <- gene}
    } else {
      orthoGroup  <- NaN
    }
    if (length(orthoGene)>0){
      orthoGene <- strsplit(orthoGene,' ')[[1]][1]
    } else{orthoGene <- NaN}
    #Search uniprot info
    uni_Index <- grepl(gene,uniprotData[,3])
    if (length(uni_Index)>0){
      uni_Index <- which(uni_Index==TRUE)[1]
      protID    <- uniprotData[uni_Index,1]
      protName  <- uniprotData[uni_Index,2]
      shortName <- uniprotData[uni_Index,4]
      MWeigth   <- uniprotData[uni_Index,6]
      SeqLength <- nchar(uniprotData[uni_Index,7])
      if (all(organism=='yli')){orthoGene<-protID}
    } else{
      protID    <- ''
      protName  <- ''
      shortName <- ''
      MWeigth   <- ''
      SeqLength <- ''
    }
    #Search gene GOterm
    GOterms  <- ''
    GO_Index <- which(GO_terms[,1] == orthoGene)
    if (length(GO_Index)>0){
      GOterms  <- GO_terms[GO_Index,2]
      GOterms  <- paste(GOterms, collapse = '//')
    }
    #Combine variables for a row in the integrated table
    #if (significance & !is.na(P_ref)){  
    if (justDEgenes) {
      conditional <- significance
      } else {conditional <- TRUE}
     
    if (conditional){  
      row       <- cbind(originalGene,condition,significance,direction,P_ref,P_cond,orthoGroup,protName,shortName,MWeigth,SeqLength,GOterms) 
      tableData <- rbind(tableData,row)
    } 
  }
}
setwd(resultsPath)
if (justDEgenes) {
  fileName <- paste(organism,'_integratedTable_DEgenes.txt',sep='')
} else {
  fileName <- paste(organism,'_integratedTable_Allgenes.txt',sep='')
}
write.table(tableData, file = fileName, row.names = F,quote = FALSE,sep="\t")
}
