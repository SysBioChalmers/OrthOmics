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
repoPath  <- '/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis'
#Internal functions
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
setwd(scriptsPath)
#Provide organism code [Sce,Kma,Yli]
organism    <- 'kma'
#dataPath    <- paste(repoPath,'/RNA-seq',sep='')
#resultsPath <- paste(dataPath,'/',organism,'/Results',sep='')
dataPath <- paste(repoPath,'/Final_products_Tolerance/Docosanol/RNAseq',sep = '')
resultsPath <- paste(repoPath,'/Final_products_Tolerance/Docosanol/Results/',organism,sep='')
#================== 1. Load data and add grouping info ====================================
setwd(scriptsPath)
source('loadRNAdata_Doc.R')
output      <- loadRNAdata_Doc(dataPath,organism)
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
output   <- filterData(dataset,replicates,'median')
filtered <- output[[1]]
detected <- output[[2]]
rm(output)
filtered.data <- dataset[filtered,]
#============ Get venn diagram for measured RNA
setwd(scriptsPath)
source('plotVennDiagram.R')
setwd(resultsPath)
png(paste(organism,'_RNA_vennAllconds.png',sep=''),width = 600, height = 600)
  intLabSize <- c(rep(2.5,3))
  intLabSize[2]<-3.5
  allConds <- plotVennDiagram(detected,conditions,colorValues,intLabSize,2)
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
output <-filterLowReads(filtered.data,lcpm2,organism,resultsPath)
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(lcpm,x,lcpm2,x2)
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
getBoxPlots(x,x2,titleStr,resultsPath,organism)
dev.off()

#================== 5. Unsupervised clustering ================================================
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)
prots.PCA <- getPCAplot(lcpm2,conditions,group,replicates,colorValues,organism)
#======================= 6. Pairwise DE analysis ==============================================
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)
x <-filtered.data
x <- DGEList(counts = x, genes = rownames(filtered.data))
#x<- data
# To properly normalize for library size we use TMM normalization, as discussed in the lectures.
x2 <- calcNormFactors(x, method = "TMM")
x2 <- estimateDisp(x2)
x2$samples$group <- group
#Call DE analysis internal function
#Define DE thresholds
logPval <- abs(log10(0.01))
log2FC  <- 0.75
output  <- DEpairwiseAnalysis(x2,organism,conditions,colorValues,logPval,log2FC)
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up   <- output[[3]]
Excsv_down <- output[[4]]