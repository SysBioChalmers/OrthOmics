RNAseqAnalysis <- function(organism,stringent,normMethod,minLCPM,logPval,log2FC,adjustedP,repoPath){
#RNAseqAnalysis
#
#Function that performs differential expression analysis on RNAseq datasets the function follows a pipeline in 
#which noisy measurements are removed, then low reads are removed and dataset is normalized. 
#PCA is also run and results stored as a plot. DE analysis is carried out using limma and edgeR.
#
# organism      (string) Organism ID (sce, kma or yli)
# stringent     TRUE if noisy measurements should be filtered out according to: datum>=(stddev(row)/median(row)). 
#               FALSE if datum>=(stddev(row)/mean(row)) should be used instead.
# normMethod    (string) recommended 'TMM'
# minLCPM       (double) minimum value for log2CPM when filtering out low reads  
# logPval       (double) abs(Log10) for the DE pValue threshold
# log2FC        (double) abs(log2FC) for the DE fold-change threshold
# adjustedP     TRUE if adjusted pValue computation should be used
# repoPath      Main repository directory
#
# Usage: RNAseqAnalysis(organism,stringent,normMethod,minLCPM,logPval,log2FC,adjustedP,repoPath)
#
# Last modified: Ivan Domenzain. 2019-05-20
#
  
#================================ RELEVANT DIRECTORIES ====================================
#Relevant paths (The user should provide the path in which the repository is stored)
#Internal functions path
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')
DBpath      <- paste(repoPath,'/Databases/Uniprot/',sep='')
#Original data path
dataPath    <- paste(repoPath,'/RNA-seq',sep='')
#Directory for results
resultsPath <- paste(dataPath,'/',organism,'/Results',sep='')
cat("#==================               1. Load data and add grouping info             ================\n")
setwd(scriptsPath)
source('loadRNAdata.R')
cat("Loading RNAseq data\n")
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
cat("#==================                      2. Filter Data                          ================\n")
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
cat(paste("Getting overlap of measured transcripts across conditions: ",dataSource,"\n",sep=""))
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
cat("#================== 3. visualize data samples distributions and filter low reads =================\n")
setwd(scriptsPath)
source('filterLowReads.R')
source('plotDistributions.R')
setwd(resultsPath)
#Plot reads dritributions (in log2CPM) for filtered and unfiltered data, the plots show a vertical 
#line in the 1 CPM threshold to indicate the region of the distributions to remove with the 
#filterLowReads script
cat(paste("Filtering low reads (log2(cpm) < 0)\n",sep=""))
x  <- dataset
x2 <- filtered.data
png(paste(organism,'_SamplesDistributions.png',sep=''),width = 1200, height = 600)
plotDistributions(x,x2,' RNA', 0.4,colorValues,replicates)
dev.off()
rm(x,x2)
#Filter low reads
indexes       <- filterLowReads(lcpm,coverage,'-',replicates,minLCPM)
filtered.data <- filtered.data[indexes,]
lcpm          <- lcpm[indexes,]
cat("#==================                   4. Data normalization                      ================\n")
setwd(scriptsPath)
source('getBoxPlots.R')
cat(paste("Data normalization using ",normMethod,' method\n',sep=""))
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
cat("#==================               5. Unsupervised clustering                     ================\n")
#Get PCA for the filtered data
setwd(scriptsPath)
source('getPCAplot.R')
setwd(resultsPath)
cat(paste("PCA for RNAseq data\n",sep=""))
plot_name <- paste(organism,'_RNAseq_PCA.png',sep='')
prots.PCA <- getPCAplot(x2$counts,conditions,grouping,replicates,colorValues,organism,plot_name,' RNA')
rm(x2)
cat("#=======================            6. Pairwise DE analysis                      ================\n")
setwd(scriptsPath)
source('DEpairwiseAnalysis.R')
setwd(resultsPath)
#Call DE analysis internal function
output   <- DEpairwiseAnalysis(filtered.data,organism,conditions,colorValues,logPval,log2FC,adjustedP,'RNA',grouping,normMethod)
upReg_AllConds   <- output[[1]]
downReg_AllConds <- output[[2]]
Excsv_Up         <- output[[3]]
Excsv_down       <- output[[4]]
cat("\n")
cat(paste("Transcriptomics analysis completed, results are available in RNA-seq/",organism,"/Results\n",sep=""))
cat("\n")
}