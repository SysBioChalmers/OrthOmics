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
organism    <- 'yli'
#Indicate if the dataset is absolute quantification of proteins
absolute <- FALSE
#================== 1. Load data and add grouping info ====================================
DBpath <- paste(repoPath,'/Databases/Uniprot/',sep='')
setwd(scriptsPath)
source('load_ProtData.R')
source('normalize_SCounts.R')
if (absolute){
  dataPath    <- paste(repoPath,'/Proteomics/Absolute/data',sep='')
  resultsPath <- paste(repoPath,'/Proteomics/Absolute/Results/',organism,sep='')
  DBpath      <- paste(repoPath,'/Databases/Uniprot/',sep='')
  #Load NSAF data
  output_1    <- load_ProtData (dataPath,c(),organism,'NSAF')
  #Load IBAQ data
  output_2    <- load_ProtData (dataPath,c(),organism,'IBAQ')
  legends <- c('NSAF','IBAQ')
} else {
  dataPath    <- paste(repoPath,'/Proteomics/Relative/data',sep='')
  resultsPath <- paste(repoPath,'/Proteomics/Relative/Results/',organism,sep='')
  #Load XIC data
  output_1      <- load_ProtData (dataPath,c(),organism,'XIC')
  #Convert to linear scale
  output_1[[1]] <- as.data.frame(10^(output_1[[1]]))
  #Load Scounts data
  output_2      <- load_ProtData (dataPath,DBpath,organism,'SCounts')
  legends <- c('XIC','SCounts')
}
dataset_1     <- output_1[[1]]
dataset_2     <- output_2[[1]]
#Get data IDs
proteins_1    <- output_1[[2]]
genes_1       <- output_1[[3]]
proteins_2    <- output_2[[2]]
genes_2       <- output_2[[3]]
#Set rownames for the dataset (proteins or genes)
rownames(dataset_1) <- genes_1
rownames(dataset_2) <- genes_2
#Get grouping information
conditions  <- output_1[[4]]
colorValues <- output_1[[5]]
replicates  <- output_1[[6]]
group       <- output_1[[7]]

rm(output_1)
rm(output_2)
rm(dataPath)
#================== 2. Filter Data ================================================
#Filter data: Keep those proteins that were measured in at least 2/3 of the replicates for at least one condition
#Remove those proteins which show a RSD>1 across triplicates for the conditions in which it was measured
#Remove those proteins with a SD == 0 across all samples
setwd(scriptsPath)
source('filterData.R')
#Coverage means the proportion of replicates in which a protein should be present 
#in order to be considered as measured for a given condition
coverage <- 2/3
#Filter by median value for XIC or NSAF
output   <- filterData(dataset_1,replicates,'median','prots',coverage)
filtered <- output[[1]]
detected_1 <- output[[2]]
rm(output)
filtered_1 <- dataset_1[filtered,]
#Filter by mean value for SC or IBAQ
output     <- filterData(dataset_2,replicates,'mean','prots',coverage)
filtered   <- output[[1]]
detected_2 <- output[[2]]
rm(output)
filtered_2 <- dataset_2[filtered,]
rm(filtered)
#============ Get venn diagram for measured genes
setwd(scriptsPath)
source('plotVennDiagram.R')
setwd(resultsPath)
#XIC
png(paste(organism,'_1_vennAllconds.png',sep=''),width = 600, height = 600)
if (all(organism == 'yli')) {
  intLabSize <- c(rep(2,7))
  intLabSize[2]<-2.5
  intLabSize[4]<-2.5
  intLabSize[6]<-2.5
  intLabSize[5]<-3
  allConds <- plotVennDiagram(detected_1,conditions,colorValues,intLabSize,3)
}else  {
  intLabSize   <- c(rep(2.5,15))
  intLabSize[6]<- 4
  allConds <- plotVennDiagram(detected_1,conditions,colorValues,intLabSize,4)
}
dev.off()
#spectral Counts
png(paste(organism,'_2_vennAllconds.png',sep=''),width = 600, height = 600)
if (all(organism == 'yli')) {
  intLabSize <- c(rep(2,7))
  intLabSize[2]<-2.5
  intLabSize[4]<-2.5
  intLabSize[6]<-2.5
  intLabSize[5]<-3
  allConds <- plotVennDiagram(detected_2,conditions,colorValues,intLabSize,3)
}else  {
  intLabSize   <- c(rep(2.5,15))
  intLabSize[6]<- 4
  allConds <- plotVennDiagram(detected_2,conditions,colorValues,intLabSize,4)
}
dev.off()
#================== 3. visualize data samples distributions and filter low reads =================
setwd(scriptsPath)
source('filterLowReads.R')
source('plotDistributions.R')
setwd(resultsPath)
#XIC
x  <- filtered_1
x <- cpm(x,log=FALSE)
#Filter low reads (log2cpm<0)
filtered_1 <-filterLowReads(filtered_1,1,'-',replicates)
x2 <- filtered_1
x2 <- cpm(x2,log=FALSE)
x  <- log10(x)
x2 <- log10(x2)
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_1_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,x2,' proteins',1)
dev.off()
rm(x,x2)

#spectral counts
x  <- filtered_2
x <- cpm(x,log=FALSE)
#Filter low reads (log2cpm<0)
filtered_2 <-filterLowReads(filtered_2,1,'-',replicates)
x2 <- filtered_2
x2 <- cpm(x2,log=FALSE)
x  <- log10(x)
x2 <- log10(x2)
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_2_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,x2,' proteins',1)
dev.off()
rm(x,x2)
#Get overlap between methods
intLabSize <- c(rep(4,3))
png(paste(organism,'_prots_methods.png',sep=''),width = 600, height = 600)
methods <- plotVennDiagram(list(rownames(filtered_1),rownames(filtered_2)),legends,c('red','blue'),intLabSize,2)
dev.off()
#================== 4. Data normalization ================================================
setwd(scriptsPath)
source('getBoxPlots.R')
#XIC
x_1 <- DGEList(counts = (filtered_1), genes = rownames(filtered_1))
x_1$samples$group <- group
x2_1 <- calcNormFactors(x_1, method = 'TMM')
plot_name <- paste(organism,'_1_normalization.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr  <- paste(organism, '_',length(filtered_1[,1]),sep='')
getBoxPlots(x_1,x2_1,titleStr,resultsPath,organism,'prots')
dev.off()
#spectral counts
x_2 <- DGEList(counts = (filtered_2), genes = rownames(filtered_2))
x_2$samples$group <- group
x2_2 <- calcNormFactors(x_2, method = 'TMM')
plot_name <- paste(organism,'_2_Box_normalization.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr  <- paste(organism, '_',length(filtered_2[,1]),sep='')
getBoxPlots(x_2,x2_2,titleStr,resultsPath,organism,'prots')
dev.off()
rm(titleStr)

#================== 5. Unsupervised clustering ================================================
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)
#The datasets are already normalized, in the sense that they take the protein length or MW into 
#account in the definition of their metrics so, there's no need to use TMM normalized data for 
#the next tasks.

#XIC or NSAF
data <- as.data.frame(x_1$counts)
plot_name <- paste(organism,'_1_PCA.png',sep='')
prots.PCA_1 <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#spectral counts or IBAQ
data <- as.data.frame(x_2$counts)
plot_name <- paste(organism,'_2_PCA.png',sep='')
prots.PCA_2 <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#======================= 6. Pairwise DE analysis ==============================================
indexes  <- which(!is.element(rownames(filtered_2),rownames(filtered_1)))
#Merge both datasets
dataset  <- rbind(filtered_1,filtered_2[indexes,]) 
dataset <- DGEList(counts = (dataset), genes = rownames(dataset))
dataset$samples$group <- group
#dataset <- calcNormFactors(dataset, method = 'TMM')
dataset <- estimateDisp(dataset)
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)

#Call DE analysis internal function
#Define DE thresholds
logPval <- abs(log10(0.01))
#A 50% percent of fold-change
log2FC  <- 1
adjusted <- TRUE
output  <- DEpairwiseAnalysis(dataset,organism,conditions,colorValues,logPval,log2FC,adjusted,'Proteins')
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up   <- output[[3]]
Excsv_down <- output[[4]]