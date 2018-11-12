install.packages('VennDiagram')
library(VennDiagram)
library(devtools)
library(ggbiplot)
library(limma)
library(edgeR)
library(tidyverse) # Collection of useful R-packages
library(RColorBrewer)
library(snowfall)

#Provide organism code [Sce,Kma,Yli]
organism <- 'Kma'
resultsPath <- '/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Proteomics/Results'
#================== 1. Data preprocessing ================================================
#Extract numeric data
data <- emPAI[,6:17]
#Convert to linear scale
data <- log10(data+1)
#Convert missing values to 0
data[is.na(data)] <-0

# Add the grouping information (and enforce the order of the samples).
if (all(organism == 'Yli')){
  conditions <- c('Ref','HiT','LpH')
  colorValues <- c("black", "red", "#009E73")
  replicates <- c(4,3,4)
  group <- factor(c(rep(conditions[1],replicates[1]),
                    rep(conditions[2],replicates[2]),
                    rep(conditions[3],replicates[3])), levels = conditions)
  geneIDs      <- emPAI$W29_NCBI
} else {
  geneIDs      <- emPAI$NCBI_CBS6556
  if (all(organism == 'Sce')){
    temp <- data
    data[,7:9]   <- temp[,10:12]
    data[,10:12] <- temp[,7:9]
    geneIDs      <- emPAI$GeneName_S288C
  }
  conditions <- c('Ref','HiT','LpH','Osm')
  colorValues <- c("black", "red", "#009E73","blue")
  replicates <- c(3,3,3,3)
  group <- factor(c(rep(conditions[1],replicates[1]),
                    rep(conditions[2],replicates[2]),
                    rep(conditions[3],replicates[3]),
                    rep(conditions[4],replicates[4])), levels = conditions)
}
#Generate colNames
cols <- c()
for (i in 1:length(conditions)) {
  for (j in 1:replicates[i]) {
    cols <- c(cols,paste(conditions[i],'_',j,sep=''))
  }
}
colnames(data) <- cols
remove(cols)

#Filter data: Keep those proteins that were measured in at least 2/3 of the replicates for all of the conditions
#Remove those proteins which show a RSD>1 across triplicates for at least one of the conditions
#Remove those proteins with a SD == 0 across triplicates for at least one one of the conditions
rowIndexes <- seq(1,length(data[,1]),1)
maxIndex <- 0
labelStr <- c()
detected <- c()
for (i in 1:length(conditions)) {
  threshold <- replicates[i]*(2/3)
  indexes <- c()
  indexesVenn<- c()
  #detected[i] <- 0
  #Replicates indexes for the i-th condition
  colIndexes <- seq(maxIndex+1,maxIndex+replicates[i],1)
  for (j in 1:length(data[,1])) {
    rowData <- as.matrix(data[j,colIndexes])
    RSD <- sd(rowData)/mean(rowData)
    if ((sum(rowData>0)>=threshold) & RSD<=1 & sd(rowData)>0) {
      indexes <- c(indexes,j)
    }
    #Save detected proteins per condition for genreating a Venn diagram
    if (sum(rowData)>0){indexesVenn <- c(indexesVenn,j)}
      
  }
  rowIndexes <- intersect(rowIndexes,indexes)
  maxIndex <- colIndexes[length(colIndexes)]
  print(length(indexes))
  labelStr <- c(labelStr, seq(1,replicates[i],1))
  detected[[i]]<-indexesVenn
}
#Get filtered values
filtered.data <- data[rowIndexes,]
geneIDs       <- geneIDs[rowIndexes]
#============ Get venn diagram for measured proteins
area   <- c()
shared <- c()
for (i in 1:3){
  if (i==3){j<-1}
  else {j<-i+1}
area[i]   <- length(detected[[i]])
shared[i] <- length(intersect(detected[[i]],detected[[j]]))
}
#Intersect 1,2,3
shared[4] <- length(intersect(intersect(detected[[1]],detected[[2]]),detected[[3]]))

if (all(organism == 'Yli')) {
  intersectLabels <- c(rep(1,7))
  intersectLabels[5]<-2
  venn.plot<-draw.triple.venn(area[1], area[2],area[3], 
                              shared[1],shared[2],shared[3],shared[4],
                              category = conditions,scaled = F,fill = colorValues, 
                              lty = c(rep("solid",3)), cex = intersectLabels, cat.cex = 2)
}else  {
  area[4] <- length(detected[[4]])
  #Intersect 2,4
  shared[5] <- length(intersect(detected[[2]],detected[[4]]))
  #Intersect 3,4
  shared[6] <- length(intersect(detected[[3]],detected[[4]]))
  #Intersect 4,1
  shared[7] <- length(intersect(detected[[1]],detected[[4]]))
  #Intersect 1,2,4
  shared[8] <- length(intersect(intersect(detected[[1]],detected[[2]]),detected[[4]]))
  #Intersect 1,3,4
  shared[9] <- length(intersect(intersect(detected[[1]],detected[[3]]),detected[[4]]))
  #Intersect 2,3,4
  shared[10] <- length(intersect(intersect(detected[[2]],detected[[3]]),detected[[4]]))
  #Intersect 1,2,3,4
  shared[11] <- length(intersect(intersect(detected[[1]],detected[[2]]),intersect(detected[[3]],detected[[4]]))) 
  intersectLabels <- c(rep(1,15))
  intersectLabels[6]<-2
   venn.plot <- draw.quad.venn(area[1], area[2],area[3],area[4], 
                               shared[1],shared[3],shared[7],shared[2],shared[5],
                               shared[6],shared[4],shared[8],shared[9],shared[10],shared[11],
                               #n12, n13,n14, n23, n24,
                               #n34, n123, n124, n134, n234, n1234, 
                               category = conditions,scaled = F,fill = colorValues, 
                               lty = c(rep("solid",4)), cex = intersectLabels, cat.cex = 2);

}
#Save plot
setwd(resultsPath)
plot_name <- paste(organism,'_vennAllconds.png')
dev.copy(png,plot_name)
dev.off()
#================== 2. visualize data samples distributions ================================================
x    <- data
lcpm <- cpm(x, log = T)

x2 <- filtered.data
lcpm2 <- cpm(x2, log = T)
# With the following code we can make a plot showing the raw log2CPM reads, the dotted
# line indicates zero logCPM (= 1 cpm)
nsamples <- ncol(x)
# Ignore warning, some samples will have identical color.
col <- brewer.pal(nsamples, "Paired") 
par(mfrow = c(1, 2))
#Plot raw data
plot(
  density(lcpm[, 1]),
  col = col[1],
  ylim = c(0, 0.8),
  las = 2,
  main = paste("A. Raw data ", length(lcpm[,1]), " proteins"),
  xlab = "Log2-cpm"
)
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[, i])
  lines(den$x, den$y, col = col[i])
}
#Plot filtered data
plot(
  density(lcpm2[, 1]),
  col = col[1],
  ylim = c(0, 0.8),
  las = 2,
  main = paste("B. Filtered data", length(lcpm2[,1]), " proteins"),
  xlab = "Log2-cpm"
)
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm2[, i])
  lines(den$x, den$y, col = col[i])
}
#Save plot
setwd(resultsPath)
plot_name <- paste(organism,'_samplesDistributions.png')
dev.copy(png,plot_name)
dev.off()
remove(x2)
remove(lcpm)
remove(lcpm2)
#================== 3. Data normalization ================================================
#=========== First normalize on rawData
x <- DGEList(counts = data, genes = rownames(data))
#x<- data
# To properly normalize for library size we use TMM normalization, as discussed in the lectures.
x2 <- calcNormFactors(x, method = "TMM")

# We then visualize this with two graphs, before and after normalization:
par(mfrow = c(1, 2))
lcpm2 <- cpm(x, log = TRUE)
boxplot(lcpm2,  las = 2, col = col, main = "")

titleStr <- paste(organism, '_RawData',length(lcpm2[,1]), ' proteins: Unnormalised')
title(main = titleStr, ylab = "Log-cpm")
x2 <- calcNormFactors(x2)
lcpm2<- cpm(x2, log = TRUE)
boxplot(lcpm2, las = 2, col = col, main = "")
titleStr <- paste(organism, '_RawData',length(lcpm2[,1]), ' proteins: Normalised')
title(main = titleStr, ylab = "Log-cpm")
#Save plot
setwd(resultsPath)
plot_name <- paste(organism,'_RawDataBoxPlots.png')
dev.copy(png,plot_name)
dev.off()
#=========== Then repeat on filtered data
x <- DGEList(counts = filtered.data, genes = rownames(filtered.data))
#x<- data
# To properly normalize for library size we use TMM normalization, as discussed in the lectures.
x2 <- calcNormFactors(x, method = "TMM")

# We then visualize this with two graphs, before and after normalization:
par(mfrow = c(1, 2))
lcpm2 <- cpm(x, log = TRUE)
boxplot(lcpm2,  las = 2, col = col, main = "")

titleStr <- paste(organism, '_FilteredData',length(lcpm2[,1]), ' proteins: Unnormalised')
title(main = titleStr, ylab = "Log-cpm")
x2 <- calcNormFactors(x2)
lcpm2<- cpm(x2, log = TRUE)
boxplot(lcpm2, las = 2, col = col, main = "")
titleStr <- paste(organism, '_FilteredData',length(lcpm2[,1]), ' proteins: Normalised')
title(main = titleStr, ylab = "Log-cpm")
#Save plot
setwd(resultsPath)
plot_name <- paste(organism,'_FilteredDataBoxPlots.png')
dev.copy(png,plot_name)
dev.off()
#================== 4. Unsupervised clustering ================================================

#PCA 
t.filtered.data<-t(filtered.data)
prots.pca <- prcomp(t.filtered.data, center = TRUE,scale. = TRUE) 
summary(prots.pca)  
#Plot PC1 and PC2 
p<-ggbiplot(prots.pca, choices = c(1,2),var.axes = FALSE,groups = group,
            varname.size = 6,labels=labelStr,labels.size = 6,ellipse=TRUE)
p <- p + theme_bw()
#Format plot
p <- p + scale_colour_manual(values=colorValues)
titleStr <- paste(organism,'/ emPAI /',length(rowIndexes), ' proteins')
p <- p + labs(title = titleStr)
p
#Save plot
setwd(resultsPath)
plot_name <- paste(organism,'_rel_emPAI_PCA.png')
ggsave(plot_name, width = 5, height=5)

#======================= 5. Pairwise DE analysis ==============================================
x <-(10^filtered.data)-1
x <-filtered.data
x <- DGEList(counts = x, genes = geneIDs)
#x<- data
# To properly normalize for library size we use TMM normalization, as discussed in the lectures.
x2 <- calcNormFactors(x, method = "TMM")
x2 <- calcNormFactors(x2)
x2 <- estimateDisp(x2)
x2$samples$group <- group
setwd(resultsPath)
#Define DE thresholds
logPval <- abs(log10(0.05))
log2FC  <- 0.5
DEgenes <- list()
FCmat <- c()
for (i in 2:length(conditions)) {
  #Differential expression analysis
  de <- exactTest(x2, pair = c(conditions[1],conditions[i])) # Reference first!
  tt <- topTags(de, n = Inf)
  # Write CSV file
  filename <- paste(organism,'_emPAI_ref_',conditions[i],'.csv')
  write.csv(tt, file = filename, row.names = T)
  # Make volcano plot
  Pval <- tt$table$PValue
  FChange <- 2^tt$table$logFC
  FCmat<- cbind(FCmat,FChange)
}
#FCmat <- as.data.frame(FCmat)  
#rownames(FCmat) <- geneIDs
#PCA 
colnames(FCmat) <- labelStr
FCmat<-t(FCmat)
labelStr <- c()
for (i in 2:length(conditions)){
  labelStr <- cbind(labelStr,paste(conditions[1],' Vs. ',conditions[i]))
}
prots.pca <- prcomp(FCmat, center = TRUE,scale. = TRUE) 
summary(prots.pca)

#======== Extract PC coefficients
library(hydroTSM)
sorted <- c()
Diff <- c()
overlap <- c()
n <- length(geneIDs)
coefficients <- prots.pca$rotation
rownames(coefficients) <- geneIDs
abs.coeffs <- abs(coefficients)
for (i in 1:length(abs.coeffs[1,])){
  column <-as.matrix(abs.coeffs[,i])
  IDs <- (geneIDs[order(-column)])
  coeffs <- (column[order(-column)])
  sorted[[i]] <- IDs[1:n]
  Diff[[i]] <- (coeffs[1:n])
}
overlap <- sorted[[1]]
#for (i in 2:length(abs.coeffs[1,])){
 overlap <- intersect(overlap,sorted[[2]])

#}
sorted <- as.data.frame(sorted)
Diff <- as.data.frame(Diff)
p<-image(as.matrix(t(Diff[,1:2])))#, ColorRamp="Days", ncolors = 70)



