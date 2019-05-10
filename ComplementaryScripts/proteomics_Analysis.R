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
organism    <- 'sce'
#Indicate the dataset that should be used foer DE analysis
selectedDataset <- 1 #1 for XIC or NSAF, 2 for Scounts or iBAQ and 3 for merged datasets
#Define DE thresholds
logPval    <- abs(log10(0.01))
log2FC     <- 0.75
adjustedP  <- TRUE
#Should the initial dataset be normalized by MW of proteins
normByMW <- TRUE
#Filter type for determination of present and non-noisy proteins in the dataset (TRUE if filter
#criteria should be applied to all conditions, FALSE if just the reference is desired to be 
#filtered)
stringent  <- TRUE
#Normalization method for DE analysis
normMethod <- 'TMM'
#===========================================================================================
#Relevant paths (The user should provide the path in which the repository is stored)
repoPath    <- '/Users/ivand/Documents/GitHub/CHASSY_multiOmics_Analysis'
#Internal functions path
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
DBpath      <- paste(repoPath,'/Databases/Uniprot/',sep='')
#================== 1. Load data and add grouping info ====================================
setwd(scriptsPath)
source('load_ProtData.R')
source('normalize_Prots_MWeight.R')
source('normalize_Prots_AALength.R')
#==================== Relative data ====================================
dataPath    <- paste(repoPath,'/Proteomics/data/relative',sep='')
resultsPath <- paste(repoPath,'/Proteomics/Results/',organism,sep='')
#Load XIC data
output_1    <- load_ProtData(dataPath,DBpath,organism,'XIC',normByMW)
#Load Scounts data
output_2    <- load_ProtData(dataPath,DBpath,organism,'SCounts',normByMW)
#Load Scounts data and transform them to normalized TPA values
#In this fuction each SCounts value is divided by the AA chain length of the 
#protein and the whole dataset is then normalized taking the MWeights into 
#account
output_TPA  <- load_ProtData(dataPath,DBpath,organism,'SCounts',FALSE,TRUE)
legends     <- c('XIC','SCounts')
#Load NSAF data [umol/g protein]
dataPath    <- paste(repoPath,'/Proteomics/data/absolute',sep='')
output_abs  <- load_ProtData(dataPath,DBpath,organism,'NSAF',FALSE)
#Sort and rename outputs
dataset_1     <- output_1[[1]]
lcpm_1        <- cpm(dataset_1, log = T)
dataset_2     <- output_2[[1]]
lcpm_2        <- cpm(dataset_2, log = T)
dataset_abs   <- output_abs[[1]]
lcpm_abs      <- cpm(dataset_abs, log = T)
dataset_TPA   <- output_TPA[[1]]
lcpm_TPA      <- cpm(dataset_TPA, log = T)
#Get data IDs
proteins_1    <- output_1[[2]]
genes_1       <- output_1[[3]]
proteins_2    <- output_2[[2]]
genes_2       <- output_2[[3]]
proteins_abs  <- output_abs[[2]]
genes_abs     <- output_abs[[3]]
proteins_TPA  <- output_TPA[[2]]
genes_TPA     <- output_TPA[[3]]
#Set rownames for the dataset (proteins or genes)
rownames(dataset_1)   <- genes_1
rownames(dataset_2)   <- genes_2
rownames(dataset_abs) <- genes_abs
rownames(dataset_TPA) <- genes_TPA
#Get grouping information
conditions  <- output_1[[4]]
colorValues <- output_1[[5]]
replicates  <- output_1[[6]]
group       <- output_1[[7]]
#Write file with normalized values in absolute data subfolder
setwd(paste(repoPath,'/Proteomics/data/absolute',sep=''))
filename <- paste(organism,'_normalized_TPA.txt',sep='')
write.table(dataset_TPA, file = filename, row.names = T,quote = F,sep='\t')

rm(output_1)
rm(output_2)
rm(output_abs)
rm(output_TPA)
rm(dataPath)
#================== 2. Filter Data ================================================
#Filter data: Keep those proteins that were measured in at least 2/3 of the replicates 
#for at least one condition.
#Remove those proteins which show a RSD>1 across triplicates for the conditions in 
#which it was measured.
#Remove those proteins with a SD == 0 across all samples
setwd(scriptsPath)
source('filterData.R')
source('getMeanAbundances.R')
#Coverage means the proportion of replicates in which a protein should be present 
#in order to be considered as measured for a given condition
coverage <- 2/3
#Filter by median value for XIC datasets (they come from exponential values they're
#likely to span several orders of magnitude
output     <- filterData(dataset_1,replicates,'median',stringent,coverage)
filtered   <- output[[1]]
detected_1 <- output[[2]]
filtered_1 <- dataset_1[filtered,]
lcpm_1     <- lcpm_1[filtered,]
rm(output)
rm(filtered)
#Filter by mean value for SC, IBAQ or NSAF
output     <- filterData(dataset_2,replicates,'mean',stringent,coverage)
filtered   <- output[[1]]
detected_2 <- output[[2]]
filtered_2 <- dataset_2[filtered,]
lcpm_2     <- lcpm_2[filtered,]
rm(output)
rm(filtered)
#Filter absolute measurements
output       <- filterData(dataset_abs,replicates,'mean',stringent,coverage)
filtered     <- output[[1]]
detected_abs <- output[[2]]
filtered_abs <- dataset_abs[filtered,]
lcpm_abs     <- lcpm_abs[filtered,]
#Filter absolute measurements
output       <- filterData(dataset_TPA,replicates,'mean',stringent,coverage)
filtered     <- output[[1]]
detected_TPA <- output[[2]]
filtered_TPA <- dataset_TPA[filtered,]
lcpm_TPA     <- lcpm_TPA[filtered,]
# Write CSV file with the filtered absolute datasets
setwd(resultsPath)
filename <- paste(organism,'_abs_NSAF_filtered.txt',sep='')
getMeanAbundances(filtered_abs,group,conditions,filename)
filename <- paste(organism,'_abs_normTPA_filtered.txt',sep='')
getMeanAbundances(filtered_TPA,group,conditions,filename)
rm(output)
rm(filtered)
#============ Get venn diagram for measured genes
setwd(scriptsPath)
source('plotVennDiagram.R')
setwd(resultsPath)
if (all(organism == 'yli')) {
    intLabSize <- c(rep(2,7))
    intLabSize[2] <-2.5
    intLabSize[4] <-2.5
    intLabSize[6] <-2.5
    intLabSize[5] <-3
    ellipses      <-3
}else  {
    intLabSize    <- c(rep(2.5,15))
    intLabSize[6] <- 4
    ellipses      <- 4
}
#XIC 
png(paste(organism,'_1_vennAllconds.png',sep=''),width = 600, height = 600)
allConds <- plotVennDiagram(detected_1,conditions,colorValues,intLabSize,ellipses)
dev.off()
#spectral Counts 
png(paste(organism,'_2_vennAllconds.png',sep=''),width = 600, height = 600)
allConds <- plotVennDiagram(detected_2,conditions,colorValues,intLabSize,ellipses)
dev.off()
#NSAF (absolute)
png(paste(organism,'_absolute_vennAllconds.png',sep=''),width = 600, height = 600)
allConds <- plotVennDiagram(detected_abs,conditions,colorValues,intLabSize,ellipses)
dev.off()
#================== 3. visualize data samples distributions and filter low reads =================
setwd(scriptsPath)
source('filterLowReads.R')
source('plotDistributions.R')
setwd(resultsPath)
#XIC 
x  <- filtered_1
#Filter low reads (log2cpm<0)
indexes    <-filterLowReads(lcpm_1,1,'-',replicates,0)
filtered_1 <- filtered_1[indexes,]
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_1_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,filtered_1,' proteins',0.4,colorValues,replicates)
dev.off()
rm(x)
#spectral counts
x  <- filtered_2
#Filter low reads (log2cpm<0)
indexes    <-filterLowReads(lcpm_2,1,'-',replicates,0)
filtered_2 <- filtered_2[indexes,]
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_2_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,filtered_2,' proteins',0.4,colorValues,replicates)
dev.off()
rm(x)

#Get overlap between methods
intLabSize <- c(rep(4,3))
png(paste(organism,'_prots_relative.png',sep=''),width = 600, height = 600)
methods <- plotVennDiagram(list(rownames(filtered_1),rownames(filtered_2)),legends,c('red','blue'),intLabSize,2)
dev.off()
#Get overlap between XIC and absolute data
intLabSize <- c(rep(4,7))
png(paste(organism,'_prots_rel_Abs.png',sep=''),width = 600, height = 600)
methods <- plotVennDiagram(list(rownames(filtered_1),rownames(filtered_2),rownames(filtered_abs)),c(legends,'NSAF'),c('red','blue','grey'),intLabSize,3)
dev.off()
#================== 4. Data normalization ================================================
setwd(scriptsPath)
source('getBoxPlots.R')
#XIC
x_1               <- DGEList(counts = (filtered_1), genes = rownames(filtered_1))
x_1$samples$group <- group
#How do the distributions look before and after normalization?
x2_1              <- calcNormFactors(x_1, method = normMethod)
plot_name         <- paste(organism,'_1_normalization.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr          <- paste(organism, '_',length(filtered_1[,1]),sep='')
getBoxPlots(x_1,x2_1,titleStr,resultsPath,organism,'prots',colorValues,replicates)
dev.off()
#spectral counts
x_2               <- DGEList(counts = (filtered_2), genes = rownames(filtered_2))
x_2$samples$group <- group
x2_2              <- calcNormFactors(x_2, method = normMethod)
plot_name         <- paste(organism,'_2_Box_normalization.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr          <- paste(organism, '_',length(filtered_2[,1]),sep='')
getBoxPlots(x_2,x2_2,titleStr,resultsPath,organism,'prots',colorValues,replicates)
dev.off()
rm(titleStr)
#Absolute data
x_abs               <- DGEList(counts = (filtered_abs), genes = rownames(filtered_abs))
x_abs$samples$group <- group
x_abs_2             <- calcNormFactors(x_abs, method = normMethod)
plot_name           <- paste(organism,'_Abs_Box_normalization.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr          <- paste(organism, '_',length(filtered_abs[,1]),sep='')
getBoxPlots(x_abs,x_abs_2,titleStr,resultsPath,organism,'prots',colorValues,replicates)
dev.off()
rm(titleStr)

#================== 5. Unsupervised clustering ================================================
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)

#NOTE: the datasets are already normalized, in the sense that they take the protein length or MW
#into account in the definition of their metrics or in the normalize_Prots_MWeight.R function, 
#there's no need to use TMM normalized data for the next tasks.

#XIC
data        <- as.data.frame(x_1$counts)
plot_name   <- paste(organism,'_1_PCA.png',sep='')
prots.PCA_1 <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#spectral counts
data        <- as.data.frame(x_2$counts)
plot_name   <- paste(organism,'_2_PCA.png',sep='')
prots.PCA_2 <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#Absolute data
data          <- as.data.frame(x_abs$counts)
plot_name     <- paste(organism,'_Absolute_PCA.png',sep='')
prots.PCA_Abs <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#======================= 6. Pairwise DE analysis ==============================================
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)

#Select the dataset that should go through the DE analysis
if (selectedDataset == 3){ #Merge both datasets (methods, prioritizing to the first one)
  indexes <- which(!is.element(rownames(filtered_2),rownames(filtered_1)))
  dataset <- rbind(filtered_1,filtered_2[indexes,]) 
} else {
  if (selectedDataset == 2){dataset  <- filtered_2}
  if (selectedDataset == 1){dataset  <- filtered_1}
}

#Call DE analysis function
output <- DEpairwiseAnalysis(dataset,organism,conditions,colorValues,logPval,log2FC,adjustedP,'Proteins',group)
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up         <- output[[3]]
Excsv_down       <- output[[4]]