#createIntegratedTable <- function(dataset,grouping,metric,stringent,coverage){
#Function that integrates the results from RNAseq DE with absolute proteomics
#measurements and genes annotation.
#
# Ivan Domenzain. created 2019-03-25
#
  
#==================================== DEFINE VARIABLES =======================================
#Provide organism code [sce,kma,yli]
organism    <- 'kma'
organism2   <- organism
if (all(organism == 'yli')){
  conditions <- c('HiT','LpH')
  colorValues <- c("red", "#009E73")
} else {
  conditions <- c('HiT','LpH','Osm')
  colorValues <- c("red", "#009E73","blue")
  if (all(organism=='cpk')){organism2 <- 'sce'}
}
#Define DE thresholds
FDR        <- 0.01
logPval    <- abs(log10(FDR))
log2FC     <- 1
adjustedP  <- TRUE
#=============================== Relevant directories ======================================
#Relevant paths (The user should provide the path in which the repository is stored)
repoPath    <- '/Users/ivand/Documents/GitHub/CHASSY_multiOmics_Analysis'
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
resultsPath <- paste(repoPath,'/Integration',sep='')
uniprotPath <- paste(repoPath,'/Databases/Uniprot',sep='')
RNApath     <- paste(repoPath,'/RNA-seq',sep='')
Protpath    <- paste(repoPath,'/Proteomics/Results/',organism2,sep='')
Orthpath    <- paste(repoPath,'/Orthologs',sep='')
GOTermsPath <- paste(repoPath,'/Databases/GO_terms',sep='')
#=============================== Load data =================================================
#============= Load absolute proteomics measurements (all conditions)
setwd(Protpath)
fileName <- paste(organism2,'_abs_NSAF_filtered.csv',sep='')
prot_Abs <- read.delim(fileName, header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
proteins <- prot_Abs[,1]
prot_Abs <- prot_Abs[,2:ncol(prot_Abs)]
#============= Load gene groups (orthology) list
setwd(Orthpath)
fileName   <- paste(organism,'_orthology_groups.txt',sep='')
geneGroups <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
#Transform groups to numeric values
if (all(organism=='kma')){
  proteins <- gsub('KMXK','KMXK_',proteins)
  geneGroups[(which(geneGroups[,2]=="C+C")),2] <- 5
  geneGroups[(which(geneGroups[,2]=="K+L")),2] <- 4
  geneGroups[(which(geneGroups[,2]=="K+S")),2] <- 3
  geneGroups[(which(geneGroups[,2]=="K+Y")),2] <- 2
  geneGroups[(which(geneGroups[,2]=="K+H")),2] <- 1
} else{
  if (all(organism=='yli')){
    geneGroups[(which(geneGroups[,2]=="W+W")),2]  <- 5
    geneGroups[(which(geneGroups[,2]=="Y+B")),2]  <- 4
    geneGroups[(which(geneGroups[,2]=="Y+CH")),2] <- 3
    geneGroups[(which(geneGroups[,2]=="Y+K")),2]  <- 2
    geneGroups[(which(geneGroups[,2]=="Y+H")),2]  <- 1
  } else {
      if (all(organism=='cpk')){
        geneGroups[(which(geneGroups[,2]=="C+C")),2]  <- 5
        geneGroups[(which(geneGroups[,2]=="C+E")),2]  <- 4
        geneGroups[(which(geneGroups[,2]=="C+K")),2]  <- 3
        geneGroups[(which(geneGroups[,2]=="C+Y")),2]  <- 2
        geneGroups[(which(geneGroups[,2]=="C+H")),2]  <- 1
      } 
  }
}
#============= Load GO terms files 
setwd(GOTermsPath)
fileName <- paste(organism2,'_GOterms.txt',sep="")
GO_terms  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")

#Loop through all conditions 
tableData <- c()
for (i in 1:length(conditions)){
  condition <- conditions[i]
  #Load RNAseq DE results for the i-th condition
  setwd(paste(RNApath,'/',organism,'/Results/DE_log2FC_',log2FC,'_FDR_',FDR,sep=''))
  fileName    <- paste(organism,'_RNA_ref_',conditions[i],'.csv',sep='')
  RNA_DE_data <- read.delim(fileName, header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
  genes       <- RNA_DE_data$genes
  for (j in 1:length(genes)){
    row  <- c()
    gene <- genes[j]
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
      group     <- geneGroups[geneIndex,2]
      #Extract ortholog ID for GO terms identification
      if (all(organism == 'cpk')){
        orthoGene <- geneGroups[geneIndex,3]
      } else {
        orthoGene <- gene
      } 
    } else {
      group  <- NaN
    }
   
    #Combine variables for a row in the integrated table
    #if (significance & !is.na(P_ref)){  
    if (significance){  
      row       <- cbind(gene,condition,significance,direction,P_ref,P_cond,group) 
      tableData <- rbind(tableData,row)
    } 
  }
}
setwd(resultsPath)
fileName    <- paste(organism,'_integratedTable.txt',sep='')
write.table(tableData, file = fileName, row.names = F,quote = FALSE,sep="\t")
#}