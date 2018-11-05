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
organism    <- 'yli'
dataPath    <- paste(repoPath,'/Proteomics/Relative/data',sep='')
resultsPath <- paste(repoPath,'/Proteomics/Relative/Results/',organism,sep='')

#================== 1. Load data and add grouping info ====================================
setwd(scriptsPath)
source('loadProtData.R')
emPAI <- TRUE
output      <- loadProtData(dataPath,organism,emPAI)
dataset     <- output[[1]]
proteins    <- output[[2]]
genes       <- output[[3]]
#Set rownames for the dataset (proteins or genes)
rownames(dataset) <- genes
conditions  <- output[[4]]
colorValues <- output[[5]]
replicates  <- output[[6]]
group       <- output[[7]]
rm(output)
rm(dataPath)
#================== 2. Filter Data ================================================
#Filter data: Keep those proteins that were measured in at least 2/3 of the replicates for at least one condition
#Remove those proteins which show a RSD>1 across triplicates for the conditions in which it was measured
#Remove those proteins with a SD == 0 across all samples
setwd(scriptsPath)
source('filterData.R')
#Transform emPAI to linear values
dataset  <- log10(dataset+1)
output   <- filterData(dataset,replicates,'mean','prots')
filtered <- output[[1]]
detected <- output[[2]]
rm(output)
filtered.data <- dataset[filtered,]
#============ Get venn diagram for measured RNA
setwd(scriptsPath)
source('plotVennDiagram.R')
setwd(resultsPath)
png(paste(organism,'_emPAI_vennAllconds.png',sep=''),width = 600, height = 600)
if (all(organism == 'yli')) {
  intLabSize <- c(rep(2,7))
  intLabSize[2]<-2.5
  intLabSize[4]<-2.5
  intLabSize[6]<-2.5
  intLabSize[5]<-3
  allConds <- plotVennDiagram(detected,conditions,colorValues,intLabSize,3)
}else  {
  intLabSize   <- c(rep(2.5,15))
  intLabSize[6]<- 4
  allConds <- plotVennDiagram(detected,conditions,colorValues,intLabSize,4)
}
dev.off()
#================== 3. visualize data samples distributions and filter low reads =================
setwd(scriptsPath)
source('filterLowReads.R')
source('plotDistributions.R')
setwd(resultsPath)
x  <- dataset
x <- cpm(x,log=TRUE)
#Filter low reads (log2cpm<0)
output <-filterLowReads(filtered.data,x)
filtered.data <- output
x2 <- filtered.data
x2 <- cpm(x2,log=TRUE)
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_emPAI_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,x2,' proteins',0.6)
dev.off()
rm(x,x2)
#================== 4. Data normalization ================================================
setwd(scriptsPath)
source('getBoxPlots.R')
x               <- DGEList(counts = (filtered.data), genes = rownames(filtered.data))
x$samples$group <- group
x2              <- calcNormFactors(x, method = 'TMM')
plot_name <- paste(organism,'_emPAI_Box_unnorm.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr  <- paste(organism, '_',length(filtered.data[,1]))
getBoxPlots(x,x2,titleStr,resultsPath,organism)
dev.off()

#================== 5. Unsupervised clustering ================================================
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)
#if (all(organism=='kma')){data  <- cpm(x2, log = T)}else{data <- filtered.data}
data <- filtered.data
plot_name <- paste(organism,'_emPAI_PCA.png')
prots.PCA <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#What's the contribution of individual genes to the PC's
par(mar=c(5.1,5.1,3.1,2.1));plot(prots.PCA$rotation,type='p',pch=19,col='black',cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,cex=1.5,
     main='PCA Loadings')

# Is there a way to interact with that graph?
library(plotly)
plot_ly(data.frame(prots.PCA$rotation),x=~PC1,y=~PC2,type='scatter',mode='markers',text=rownames(filtered.data),
        marker = list(size = 12,color = 'rgba(0,0,0,1)'))%>%
  layout(title = 'PCA Loadings',
         titlefont=list(family="arial",size=24,color='black'),
         font=list(family = "arial",size = 18,color = 'black'),
         margin=list(l=55,r=50,t=70,b=53,pad=0),
         hoverlabel=list(bgcolor='white',bordercolor='black',font=list(family="arial",size=16,color='black'))
  )
#======================= 6. Pairwise DE analysis ==============================================
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)
x2 <- estimateDisp(x2)
#Call DE analysis internal function
#Define DE thresholds
logPval <- abs(log10(1))
#A 50% percent of fold-change
log2FC  <- 0.5
adjusted <- FALSE
output  <- DEpairwiseAnalysis(x2,organism,conditions,colorValues,logPval,log2FC,adjusted,'Proteins')
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up   <- output[[3]]
Excsv_down <- output[[4]]