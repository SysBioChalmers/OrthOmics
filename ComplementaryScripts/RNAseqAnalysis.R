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
repoPath  <- '/Users/ivand/Documents/GitHub/CHASSY_multiOmics_Analysis'
#Internal functions
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
setwd(scriptsPath)
#Provide organism code [Sce,Kma,Yli]
organism    <- 'sce'
dataPath    <- paste(repoPath,'/RNA-seq',sep='')
resultsPath <- paste(dataPath,'/',organism,'/Results',sep='')

#================== 1. Load data and add grouping info ====================================
setwd(scriptsPath)
source('loadRNAdata.R')
output      <- loadRNAdata(dataPath,organism)
dataset     <- output[[1]]
conditions  <- output[[2]]
colorValues <- output[[3]]
replicates  <- output[[4]]
group       <- output[[5]]
rm(output)
rm(dataPath)
#================== 2. Filter Data ================================================
#Filter data: Keep those RNA that were measured in at least 2/3 of the replicates for at least one condition
#Remove those RNA which show a RSD>1 across triplicates for the conditions in which it was measured
#Remove those RNA with a SD == 0 across all samples
setwd(scriptsPath)
source('filterData.R')
output   <- filterData(dataset,replicates,'mean','RNA')
filtered <- output[[1]]
detected <- output[[2]]
rm(output)
filtered.data <- dataset[filtered,]
#============ Get venn diagram for measured RNA
setwd(scriptsPath)
source('plotVennDiagram.R')
setwd(resultsPath)
png(paste(organism,'_RNA_vennAllconds.png',sep=''),width = 600, height = 600)
if (all(organism == 'yli')) {
  intLabSize <- c(rep(3,7))
  intLabSize[5]<-4
  allConds <- plotVennDiagram(detected,conditions,colorValues,intLabSize,3)
}else  {
  intLabSize   <- c(rep(3,15))
  intLabSize[6]<- 4
  allConds <- plotVennDiagram(detected,conditions,colorValues,intLabSize,4)
}
dev.off()
#================== 3. visualize data samples distributions and filter low reads =================
setwd(scriptsPath)
source('filterLowReads.R')
source('plotDistributions.R')
setwd(resultsPath)
x      <- dataset
lcpm   <- cpm(x, log = T)
x2     <- filtered.data
lcpm2  <- cpm(x2, log = T)
#Filter low erads (log2cpm<0)
output <-filterLowReads(filtered.data,lcpm2,'-')
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(lcpm,lcpm2,' RNA', 0.3)
dev.off()
filtered.data <- output
rm(lcpm,x,x2)
#================== 4. Data normalization ================================================
setwd(scriptsPath)
source('getBoxPlots.R')
x               <- DGEList(counts = (filtered.data), genes = rownames(filtered.data))
x$samples$group <- group
x2              <- calcNormFactors(x, method = "TMM")
plot_name <- paste(organism,'_RNA_Box_unnorm.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr  <- paste(organism, '_',length(filtered.data[,1]), ' RNA: Unnormalised')
getBoxPlots(x,x2,titleStr,resultsPath,organism,'RNA')
dev.off()

#================== 5. Unsupervised clustering ================================================
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)
if (all(organism=='kma')){ 
  data  <- cpm(x2, log = T)
}else{
  data <- filtered.data # I think all organisms should use the cpm trasnformation of x2, since this is TMM normalized while filtered.data is not. This may require x3 <- cpm(filtered.data, normalized.lib.sizes = TRUE) and as.data.frame(x3)
  }
plot_name <- paste(organism,'_RNAseq_PCA.png',sep='')
prots.PCA <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' RNA')
#======================= 6. Pairwise DE analysis ==============================================
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)

x2 <- estimateDisp(x2)
#Call DE analysis internal function
#Define DE thresholds
logPval  <- abs(log10(0.01))
log2FC   <- 0.5
adjusted <- TRUE
output   <- DEpairwiseAnalysis(x2,organism,conditions,colorValues,logPval,log2FC,adjusted,'RNA')
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up   <- output[[3]]
Excsv_down <- output[[4]]