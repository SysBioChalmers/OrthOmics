proteomics_Analysis <- function(organism,dataSource,normByMW,stringent,normMethod,logPval,log2FC,adjustedP,repoPath){
#proteomics_Analysis
#
#Function that performs differential expression analysis on proteomics datasets (absolute and relative)
#the function follows a pipeline in which noisy measurements are removed, then low reads are removed and dataset 
#is normalized. PCA is also run and results stored as a plot. DE analysis is carried out using limma and edgeR.
#
# organism      (string) Organism ID (sce, kma or yli)
# dataSource    (string) Proteomics measurement dataset that should be used for DE analysis
#               (XIC, SCounts, NSAF or iBAQ)
# normByMW      TRUE if dataset should be normalized taking MWs into account
# stringent     TRUE if noisy measurements should be filtered out according to: datum>=(stddev(row)/median(row)). 
#               FALSE if datum>=(stddev(row)/mean(row)) should be used instead.
# normMethod    (string) recommended 'TMM'
# logPval       (double) abs(Log10) for the DE pValue threshold
# log2FC        (double) abs(log2FC) for the DE fold-change threshold
# adjustedP     TRUE if adjusted pValue computation should be used
# repoPath      Main repository directory
#
# Usage: proteomics_Analysis(organism,dataSource,normByMW,stringent,normMethod,logPval,log2FC,adjustedP,repoPath)
#
# Last modified: Ivan Domenzain. 2019-05-20
#
  
#Internal functions path
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
DBpath      <- paste(repoPath,'/Databases/Uniprot/',sep='')
cat("#==================               1. Load data and add grouping info             ================\n")
setwd(scriptsPath)
source('load_ProtData.R')
source('normalize_Prots_MWeight.R')
source('normalize_Prots_AALength.R')
source('sumProteinMass.R')
#==================== Relative data ====================================
dataPath    <- paste(repoPath,'/Proteomics/data/relative',sep='')
resultsPath <- paste(repoPath,'/Proteomics/Results/',organism,sep='')
cat("Loading relative proteomics data\n")
#Load Relative data
output_1    <- load_ProtData(dataPath,DBpath,organism,dataSource,normByMW)
#Load NSAF data [umol/g protein]
cat("Loading semi-absolute proteomics data [NSAF]\n")
dataPath    <- paste(repoPath,'/Proteomics/data/absolute',sep='')
output_abs  <- load_ProtData(dataPath,DBpath,organism,'NSAF',FALSE)
#Sort and rename outputs
dataset_1     <- output_1[[1]]
lcpm_1        <- cpm(dataset_1, log = T)
dataset_abs   <- output_abs[[1]]
lcpm_abs      <- cpm(dataset_abs, log = T)
#Get data IDs
proteins_1    <- output_1[[2]]
genes_1       <- output_1[[3]]
proteins_abs  <- output_abs[[2]]
genes_abs     <- output_abs[[3]]
#Set rownames for the dataset (proteins or genes)
rownames(dataset_1)   <- genes_1
rownames(dataset_abs) <- genes_abs
#Get grouping information
conditions  <- output_1[[4]]
colorValues <- output_1[[5]]
replicates  <- output_1[[6]]
group       <- output_1[[7]]
rm(output_1)
rm(output_abs)
rm(dataPath)

cat("#==================                      2. Filter Data                          ================\n")
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
cat(paste("Filtering ",dataSource," data: removal of noisy measurements\n",sep=""))
output     <- filterData(dataset_1,replicates,'median',stringent,coverage)
filtered   <- output[[1]]
detected_1 <- output[[2]]
filtered_1 <- dataset_1[filtered,]
lcpm_1     <- lcpm_1[filtered,]
rm(output)
rm(filtered)
#Filter absolute measurements
cat("Filtering NSAF data: removal of noisy measurements\n")
output       <- filterData(dataset_abs,replicates,'mean',FALSE,coverage)
filtered     <- output[[1]]
detected_abs <- output[[2]]
filtered_abs <- dataset_abs[filtered,]
proteins_abs <- proteins_abs[filtered]
lcpm_abs     <- lcpm_abs[filtered,]
# Write CSV file with the filtered absolute datasets
setwd(resultsPath)
filename      <- paste(organism,'_abs_NSAF_filtered.txt',sep='')
meanABSvalues <- getMeanAbundances(filtered_abs,group,conditions,filename)
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
#Relative 
cat(paste("Getting overlap of measured proteins across conditions: ",dataSource,"\n",sep=""))
png(paste(organism,'_',dataSource,'_vennAllconds.png',sep=''),width = 600, height = 600)
allConds <- plotVennDiagram(detected_1,conditions,colorValues,intLabSize,ellipses)
dev.off()
#NSAF (absolute)
#cat("Getting overlap of measured proteins across conditions: NSAF\n")
#png(paste(organism,'_NSAF_vennAllconds.png',sep=''),width = 600, height = 600)
#allConds <- plotVennDiagram(detected_abs,conditions,colorValues,intLabSize,ellipses)
#dev.off()

cat("#================== 3. visualize data samples distributions and filter low reads =================\n")
setwd(scriptsPath)
source('filterLowReads.R')
source('plotDistributions.R')
setwd(resultsPath)
#Relative 
cat(paste("Filtering low reads (log2(cpm) < 0) for ",dataSource,"\n",sep=""))
x  <- filtered_1
#Filter low reads (log2cpm<0)
indexes    <-filterLowReads(lcpm_1,1,'-',replicates,0)
filtered_1 <- filtered_1[indexes,]
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_',dataSource,'_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,filtered_1,' proteins',0.4,colorValues,replicates)
dev.off()
rm(x)
#Absolute 
cat(paste("Filtering low reads (log2(cpm) < 0) for NSAF\n",sep=""))
x  <- filtered_abs
#Filter low reads (log2cpm<0)
indexes    <-filterLowReads(lcpm_abs,1,'-',replicates,0)
filtered_abs <- filtered_abs[indexes,]
#Plot reads dritributions for filtered and unfiltered data
png(paste(organism,'_NSAF_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,filtered_abs,' proteins',0.4,colorValues,replicates)
dev.off()
rm(x)
#Get overlap between relative and absolute data
intLabSize <- c(rep(4,3))
png(paste(organism,'_prots_rel_Abs.png',sep=''),width = 600, height = 600)
methods <- plotVennDiagram(list(rownames(filtered_1),rownames(filtered_abs)),c('Relative','Absolute'),c('red','blue'),intLabSize,2)
dev.off()

cat("#==================                   4. Data normalization                      ================\n")
setwd(scriptsPath)
source('getBoxPlots.R')
#relative
cat(paste("Data normalization using ",normMethod,' method\n',sep=""))
x_1               <- DGEList(counts = (filtered_1), genes = rownames(filtered_1))
x_1$samples$group <- group
#How do the distributions look before and after normalization?
x2_1              <- calcNormFactors(x_1, method = normMethod)
plot_name         <- paste(organism,'_',dataSource,'_normalization.png',sep='')
setwd(resultsPath)
png(plot_name,width = 900, height = 600)
titleStr          <- paste(organism, '_',length(filtered_1[,1]),sep='')
getBoxPlots(x_1,x2_1,titleStr,resultsPath,organism,'prots',colorValues,replicates)
dev.off()
x_abs               <- DGEList(counts = (filtered_abs), genes = rownames(filtered_abs))
x_abs$samples$group <- group
cat("#==================               5. Unsupervised clustering                     ================\n")
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)
#NOTE: the datasets are already normalized, in the sense that they take the protein length or MW
#into account in the definition of their metrics or in the normalize_Prots_MWeight.R function, 
#there's no need to use TMM normalized data for the next tasks.

#Relative
cat(paste("PCA for ",dataSource," data\n",sep=""))
data        <- as.data.frame(x_1$counts)
plot_name   <- paste(organism,'_',dataSource,'_PCA.png',sep='')
prots.PCA_1 <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')
#Absolute data
cat("PCA for NSAF data\n")
data          <- as.data.frame(x_abs$counts)
plot_name     <- paste(organism,'_NSAF_PCA.png',sep='')
prots.PCA_Abs <- getPCAplot(data,conditions,group,replicates,colorValues,organism,plot_name,' Proteins')

cat("#=======================            6. Pairwise DE analysis                      ================\n")
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)
#Select the dataset that should go through the DE analysis
dataset  <- filtered_1
#Call DE analysis function
output <- DEpairwiseAnalysis(dataset,organism,conditions,colorValues,logPval,log2FC,adjustedP,'Proteins',group)
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up         <- output[[3]]
Excsv_down       <- output[[4]]
cat("\n")
cat(paste("Proteomics analysis completed, results are available in proteomics/Results/", organism,"\n",sep=""))
cat("\n")
}