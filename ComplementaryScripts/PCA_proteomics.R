library(devtools)
library(ggbiplot)
organism <- 'K. marxianus'
#Extract numeric data
data <- emPAI[,6:17]
#Convert to linear scale
data <- log10(data+1)
data[is.na(data)] <-0
# Add the grouping information (and enforce the order of the samples).
conditions <- c('Ref','HiT','LpH','Osm')
replicates <- c(3,3,3,3)
temp <- data
#data[,7:9] <- temp[,10:12]
#data[,10:12] <- temp[,7:9]
group <- factor(c(rep(conditions[1],replicates[1]),
                  rep(conditions[2],replicates[2]),
                  rep(conditions[3],replicates[3]),
                  rep(conditions[4],replicates[4])), levels = conditions)
              
#Filter data
rowIndexes <- seq(1,length(data[,1]),1)
maxIndex <- 0
labelStr <- c()
for (i in 1:length(conditions)) {
  threshold <- replicates[i]*(2/3)
  indexes <- c()
  #Replicates indexes
  colIndexes <- seq(maxIndex+1,maxIndex+replicates[i],1)
  for (j in 1:length(data[,1])) {
    rowData <- as.matrix(data[j,colIndexes])
    RSD <- sd(rowData)/mean(rowData)
    if ((sum(rowData>0)>=threshold) & RSD<=1 & sd(rowData)>0) {
      indexes <- c(indexes,j)
    }
  }
  rowIndexes <- intersect(rowIndexes,indexes)
  maxIndex <- colIndexes[length(colIndexes)]
  print(length(indexes))
  labelStr <- c(labelStr, seq(1,replicates[i],1))
}
#Get filtered values
filtered.data <- data[rowIndexes,]


#================== 2 Filtering low reads ================================================

# In x$samples$lib.size we already see that the number of reads varies a bit.
# Let's normalize for library size, by taking the log count-per-million. (In lecture is
# discussed how (log)CPM is not a satisfactory normalization, but here we just apply it
# to filter low reads!).
cpm <- cpm(filtered.data)
lcpm <- cpm(filtered.data, log = T)

# With the following code we can make a plot showing the raw logCPM reads, the dotted
# line indicates zero logCPM (= 1 cpm)
nsamples <- ncol(filtered.data)
col <-
  brewer.pal(nsamples, "Paired") # Ignore warning, some samples will have identical color.
par(mfrow = c(1, 2))
plot(
  density(lcpm[, 1]),
  col = col[1],
  ylim = c(0, 0.25),
  las = 2,
  main = "A. Raw data",
  xlab = "Log-cpm"
)
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[, i])
  lines(den$x, den$y, col = col[i])
}




#PCA 
filtered.data<-t(filtered.data)
prots.pca <- prcomp(filtered.data, center = TRUE,scale. = TRUE) 
summary(prots.pca)  
#Plot PC1 and PC2 
p<-ggbiplot(prots.pca, choices = c(1,2),var.axes = FALSE,groups = group,
            varname.size = 6,labels=labelStr,labels.size = 6,ellipse=TRUE)
p <- p + theme_bw()
#Format plot
p <- p + scale_colour_manual(values=c("black", "red", "#009E73","blue"))
#p <- p + scale_colour_discrete(name  ="Conditions")
titleStr <- paste(organism,'/ emPAI /',length(rowIndexes), ' proteins')
p <- p + labs(title = titleStr)
p
#Save plot
path <- paste('/Users/ivand/Box/CHASSY/DataAnalysis/Proteomics/figures/AllRelative_DE')
setwd(path)
plot_name <- paste(organism,'_rel_emPAI_PCA.png')
ggsave(plot_name, width = 5, height=5)