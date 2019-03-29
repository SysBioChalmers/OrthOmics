install.packages('VennDiagram')
library(VennDiagram)
library(devtools)
library(ggplot2)
library(limma)
library(edgeR)
library(tidyverse) # Collection of useful R-packages
library(RColorBrewer)
library(snowfall)
library(ggbiplot)
library(ggplot2)
#====================================DEFINE VARIABLES =======================================
#Provide organism code [sce,kma,yli]
organism        <- 'kma'
#Filter type for determination of present and non-noisy proteins in the dataset (TRUE if filter
#criteria should be applied to all conditions, FALSE if just the reference is desired to be 
#filtered)
stringent       <- TRUE
#What is a low read in the context of this dataset?
minimumlogCPM   <- 0
#Normalization method for DE analysis
normMethod      <- 'TMM'
#Define DE thresholds
logPval         <- abs(log10(0.01))
log2FC          <- 0.75
adjustedP       <- TRUE
#================================ RELEVANT DIRECTORIES ====================================
#Relevant paths (The user should provide the path in which the repository is stored)
repoPath    <- '/Users/ivand/Documents/GitHub/CHASSY_multiOmics_Analysis'
#Internal functions path
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
DBpath      <- paste(repoPath,'/Databases/Uniprot/',sep='')
#Original data path
dataPath    <- paste(repoPath,'/RNA-seq',sep='')
#Directory for results
resultsPath <- paste(dataPath,'/',organism,'/Results',sep='')
#================== 1. Load data and add grouping info ====================================
setwd(scriptsPath)
source('loadRNAdata.R')
output            <- loadRNAdata(dataPath,organism)
dataset           <- output[[1]]
lcpm              <- cpm(dataset, log = T)
genes             <- output[[2]]
conditions        <- output[[3]]
colorValues       <- output[[4]]
replicates        <- output[[5]]
grouping          <- output[[6]]
rownames(dataset) <- genes
rm(output)
rm(dataPath)
#================== 2. Filter Data ================================================
#Filter data: Keep those RNA that were measured in at least (coverage) of the replicates 
#for at least one condition remove those RNA which show a RSD>1 across triplicates for 
#the conditions in which it was measured remove those RNA with a SD == 0 across all 
#samples
setwd(scriptsPath)
source('filterData.R')
#Coverage means the proportion of replicates in which a transcript should be present 
#in order to be considered as measured for a given condition
coverage      <- 2/3
output        <- filterData(dataset,replicates,'mean',stringent,coverage)
filtered      <- output[[1]]
detected      <- output[[2]]
rm(output)
filtered.data <- dataset[filtered,]
lcpm          <- lcpm[filtered,]
#============ Get venn diagram for measured RNA
setwd(scriptsPath)
source('plotVennDiagram.R')
setwd(resultsPath)
png(paste(organism,'_RNA_vennAllconds.png',sep=''),width = 600, height = 600)
if (all(organism == 'yli')) {
  intLabSize   <- c(rep(3,7))
  intLabSize[5]<- 4
  ellipses     <- 3
  allConds     <- plotVennDiagram(detected,conditions,colorValues,intLabSize,ellipses)
}else  {
  intLabSize   <- c(rep(3,15))
  intLabSize[6]<- 4
  ellipses     <- 4
  allConds     <- plotVennDiagram(detected,conditions,colorValues,intLabSize,ellipses)
}
dev.off()
#================== 3. visualize data samples distributions and filter low reads =================
setwd(scriptsPath)
source('filterLowReads.R')
source('plotDistributions.R')
setwd(resultsPath)
#Plot reads dritributions (in log2CPM) for filtered and unfiltered data, the plots show a vertical 
#line in the 1 CPM threshold to indicate the region of the distributions to remove with the 
#filterLowReads script
x  <- dataset
x2 <- filtered.data
png(paste(organism,'_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,x2,' RNA', 0.3,colorValues,replicates)
dev.off()
rm(x,x2)
#Filter low reads
indexes       <- filterLowReads(lcpm,coverage,'-',replicates,minimumlogCPM)
filtered.data <- filtered.data[indexes,]
lcpm          <- lcpm[indexes,]
#================== 4. Data normalization ================================================
setwd(scriptsPath)
source('getBoxPlots.R')
x               <- DGEList(counts = (filtered.data), genes = rownames(filtered.data))
x$samples$group <- grouping
x2              <- calcNormFactors(x, method = "TMM")
plot_name <- paste(organism,'_RNA_Box_normalization.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr  <- paste(organism, '_',length(filtered.data[,1]), ' RNA: Unnormalised')
getBoxPlots(x,x2,titleStr,resultsPath,organism,'RNA',colorValues,replicates)
dev.off()
rm(x)
#================== 5. Unsupervised clustering ================================================
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)
plot_name <- paste(organism,'_RNAseq_PCA.png',sep='')
prots.PCA <- getPCAplot(x2$counts,conditions,grouping,replicates,colorValues,organism,plot_name,' RNA')
rm(x2)
#======================= 6. Pairwise DE analysis ==============================================
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)
#Call DE analysis internal function
output   <- DEpairwiseAnalysis(filtered.data,organism,conditions,colorValues,logPval,log2FC,adjustedP,'RNA',grouping,normMethod)
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up         <- output[[3]]
Excsv_down       <- output[[4]]