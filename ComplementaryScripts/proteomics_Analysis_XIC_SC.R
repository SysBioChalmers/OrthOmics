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
organism    <- 'kma'
dataPath    <- paste(repoPath,'/Proteomics/Relative/data',sep='')
resultsPath <- paste(repoPath,'/Proteomics/Relative/Results/',organism,sep='')

#================== 1. Load data and add grouping info ====================================
setwd(scriptsPath)
source('load_XICData.R')
#Load XIC data
output      <- load_XICData (dataPath,organism,'XIC')
dataset_XIC     <- as.data.frame(10^(output[[1]]))
proteins_XIC    <- output[[2]]
genes_XIC       <- output[[3]]
#Set rownames for the dataset (proteins or genes)
rownames(dataset_XIC) <- genes_XIC
conditions  <- output[[4]]
colorValues <- output[[5]]
replicates  <- output[[6]]
group       <- output[[7]]
rm(output)
#Load Scounts data
output      <- load_XICData (dataPath,organism,'SCounts')
dataset_SC     <- output[[1]]
proteins_SC    <- output[[2]]
genes_SC       <- output[[3]]
#Set rownames for the dataset (proteins or genes)
rownames(dataset_SC) <- genes_SC
rm(output)
rm(dataPath)
#================== 2. Filter Data ================================================
#Filter data: Keep those proteins that were measured in at least 2/3 of the replicates for at least one condition
#Remove those proteins which show a RSD>1 across triplicates for the conditions in which it was measured
#Remove those proteins with a SD == 0 across all samples
setwd(scriptsPath)
source('filterData.R')
#Transform emPAI to linear values
#dataset  <- log10(dataset+1)

#Filter by median value for XIC
output   <- filterData(dataset_XIC,replicates,'median','prots')
filtered <- output[[1]]
detected_XIC <- output[[2]]
rm(output)
filtered_XIC <- dataset_XIC[filtered,]
#Filter by mean value for SC
output   <- filterData(dataset_SC,replicates,'mean','prots')
filtered <- output[[1]]
detected_SC <- output[[2]]
rm(output)
filtered_SC <- dataset_SC[filtered,]
rm(filtered)
#============ Get venn diagram for measured RNA
setwd(scriptsPath)
source('plotVennDiagram.R')
setwd(resultsPath)
#XIC
png(paste(organism,'_XIC_vennAllconds.png',sep=''),width = 600, height = 600)
if (all(organism == 'yli')) {
  intLabSize <- c(rep(2,7))
  intLabSize[2]<-2.5
  intLabSize[4]<-2.5
  intLabSize[6]<-2.5
  intLabSize[5]<-3
  allConds <- plotVennDiagram(detected_XIC,conditions,colorValues,intLabSize,3)
}else  {
  intLabSize   <- c(rep(2.5,15))
  intLabSize[6]<- 4
  allConds <- plotVennDiagram(detected_XIC,conditions,colorValues,intLabSize,4)
}
dev.off()
#spectral Counts
png(paste(organism,'_SCounts_vennAllconds.png',sep=''),width = 600, height = 600)
if (all(organism == 'yli')) {
  intLabSize <- c(rep(2,7))
  intLabSize[2]<-2.5
  intLabSize[4]<-2.5
  intLabSize[6]<-2.5
  intLabSize[5]<-3
  allConds <- plotVennDiagram(detected_SC,conditions,colorValues,intLabSize,3)
}else  {
  intLabSize   <- c(rep(2.5,15))
  intLabSize[6]<- 4
  allConds <- plotVennDiagram(detected_SC,conditions,colorValues,intLabSize,4)
}
dev.off()
#================== 3. visualize data samples distributions and filter low reads =================
setwd(scriptsPath)
source('filterLowReads.R')
source('plotDistributions.R')
setwd(resultsPath)
#XIC
x  <- dataset_XIC
#x <- cpm(x,log=TRUE)
x <- cpm(x,log=FALSE)
x <- log10(x)
#Filter low reads (log2cpm<0)
output <-filterLowReads(filtered_XIC,x,'XIC')
filtered_XIC <- output
x2 <- filtered_XIC
x2 <- cpm(x2,log=FALSE)
x2 <- log10(x2)
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_XIC_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,x2,' proteins',1)
dev.off()
rm(x,x2)
#spectral counts
x  <- dataset_SC
x <- cpm(x,log=TRUE)
#Filter low reads (log2cpm<0)
output <-filterLowReads(filtered_SC,x,'SCounts')
filtered_SC <- output
x2 <- filtered_SC
x2 <- cpm(x2,log=TRUE)
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_SCounts_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,x2,' proteins',0.4)
dev.off()
rm(x,x2)

#================== 4. Data normalization ================================================
setwd(scriptsPath)
source('getBoxPlots.R')
#XIC
x_XIC <- DGEList(counts = (filtered_XIC), genes = rownames(filtered_XIC))
x_XIC$samples$group <- group
x2_XIC <- calcNormFactors(x_XIC, method = 'TMM')
plot_name <- paste(organism,'_XIC_Box_unnorm.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr  <- paste(organism, '_',length(filtered_XIC[,1]),sep='')
getBoxPlots(x_XIC,x2_XIC,titleStr,resultsPath,organism,'prots')
dev.off()
#spectral counts
x_SC <- DGEList(counts = (filtered_SC), genes = rownames(filtered_SC))
x_SC$samples$group <- group
x2_SC <- calcNormFactors(x_SC, method = 'TMM')
plot_name <- paste(organism,'_SCounts_Box_unnorm.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr  <- paste(organism, '_',length(filtered_SC[,1]),sep='')
getBoxPlots(x_SC,x2_SC,titleStr,resultsPath,organism,'prots')
dev.off()
rm(titleStr)

#================== 5. Unsupervised clustering ================================================
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)
#if (all(organism=='kma')){data  <- cpm(x2, log = T)}else{data <- filtered.data}
#XIC
data <- filtered_XIC
plot_name <- paste(organism,'_XIC_PCA.png',sep='')
prots.PCA <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#What's the contribution of individual genes to the PC's
par(mar=c(5.1,5.1,3.1,2.1));plot(prots.PCA$rotation,type='p',pch=19,col='black',cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,cex=1.5,
     main='PCA Loadings')
#spectral counts
data <- filtered_SC
plot_name <- paste(organism,'_SCounts_PCA.png',sep='')
prots.PCA <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#What's the contribution of individual genes to the PC's
#par(mar=c(5.1,5.1,3.1,2.1));plot(prots.PCA$rotation,type='p',pch=19,col='black',cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,cex=1.5,
#                                 main='PCA Loadings')





# Is there a way to interact with that graph?
#library(plotly)
#plot_ly(data.frame(prots.PCA$rotation),x=~PC1,y=~PC2,type='scatter',mode='markers',text=rownames(filtered.data),
#        marker = list(size = 12,color = 'rgba(0,0,0,1)'))%>%
#  layout(title = 'PCA Loadings',
#         titlefont=list(family="arial",size=24,color='black'),
#         font=list(family = "arial",size = 18,color = 'black'),
#         margin=list(l=55,r=50,t=70,b=53,pad=0),
#         hoverlabel=list(bgcolor='white',bordercolor='black',font=list(family="arial",size=16,color='black'))
#  )
#======================= 6. Pairwise DE analysis ==============================================

indexes  <- which(!is.element(rownames(filtered_SC),rownames(filtered_XIC)))
#Merge XIC and Scounts data
dataset  <- rbind(filtered_XIC,filtered_SC[indexes,]) 
dataset <- DGEList(counts = (dataset), genes = rownames(dataset))
dataset$samples$group <- group
dataset <- calcNormFactors(dataset, method = 'TMM')
dataset <- estimateDisp(dataset)
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)

#Call DE analysis internal function
#Define DE thresholds
logPval <- abs(log10(0.05))
#A 50% percent of fold-change
log2FC  <- 0.5
adjusted <- FALSE
output  <- DEpairwiseAnalysis(dataset,organism,conditions,colorValues,logPval,log2FC,adjusted,'Proteins')
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up   <- output[[3]]
Excsv_down <- output[[4]]