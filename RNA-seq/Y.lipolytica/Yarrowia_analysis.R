# R-script for Y.lipolytica dataset CHASSY

# Set working directory to data folder (note the use of / instead of \ !!)
setwd('C:/YarrowiaSeq')
# setwd('~/Documents/CHASSY')
# setwd('~/Work/CHASSY')

# Confirm we're in the write directory and look at its content:
getwd()
dir()

library(limma)
library(edgeR)
library(tidyverse) # Collection of useful R-packages
library(RColorBrewer)
library(biomaRt)
library(Glimma)
library(piano)
library(snowfall)


#================== 1. Loading and restructure data ================================================

# Import the raw counts as provided by featureCounts.
# Skip first line, contains parameters from featureCounts
counts <- read.delim('YL_STAR_feat2.txt', skip = 1)
# Read the documentation if you want to know more about a function
?read.delim
# Have a look at the table in R Studio
View(counts)

# Use gene names as row names
rownames(counts) <- counts$Geneid
# Remove other columns not containing count data
counts <- counts[,-1:-6]
# Rename columns (remove phrases first line removes beginning of title, second removes end)
colnames(counts) <- gsub('Users.doughty.Desktop.Omics_Paper.YarrowiaSeq.Chalmers_Sysbio.69875808.', '', colnames(counts))
colnames(counts) <- gsub('_STARAligned.out.bam', '', colnames(counts))


# Move reference to first columns
counts <- counts[,c(5:8,9:11,1:4)]

counts

#Make a datafile with counts if necessary
# write.csv(counts, file = "YLI_RNA_RAW_counts.csv", row.names = T)

# Move data into DGEList, which can store additional metadata
x <- DGEList(counts = counts, genes = rownames(counts))


# Add the grouping information (and enforce the order of the samples).
group <- factor(c(rep('ref',4),
                  rep('temp',3),
                  rep('pH',4)), levels = c('ref','temp','pH'))
x$samples$group <- group

# Now inspect what the data looks like:
x

# The content can also be summarized with:
dim(x)
# Showing the number of genes and the number of samples

#================== 2 Filtering low reads ================================================

# In x$samples$lib.size we already see that the number of reads varies a bit.
# Let's normalize for library size, by taking the log count-per-million. (In lecture is
# discussed how (log)CPM is not a satisfactory normalization, but here we just apply it
# to filter low reads!).
cpm <- cpm(x)
lcpm <- cpm(x, log = T)
nsamples <- ncol(x)

# With the following code we can make a plot showing the raw logCPM reads, the dotted
# line indicates zero logCPM (= 1 cpm)

col <-
  brewer.pal(nsamples, "Paired") # Ignore warning, some samples will have identical color.
par(mfrow = c(1, 2))
plot(
  density(lcpm[, 1]),
  col = col[1],
  ylim = c(0, 0.21),
  las = 2,
  main = "A. Raw data",
  xlab = "Log-cpm"
)
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[, i])
  lines(den$x, den$y, col = col[i])
}

# There are many genes with low reads. Different filters can be used to
# remove these low counts. Here we apply a filter that states that for each
# gene at least 3 of the samples should have a CPM value higher than 1. For
# the smallest library (ana.5_S19) that means 11 reads.

keep.exprs <- rowSums(cpm > 1) >= 3
x <- x[keep.exprs,, keep.lib.sizes = FALSE]
dim(x)
# This 3 means we need 3 samples or more with 1 count per million or the gene is taken from the dataset
# Now we have less genes (compared to when we ran dim(x) above)

# I DID NOT DO THIS FOR Y.LIPO 
# Alternatively, only keep those genes for which we have data in all samples:
x <- DGEList(counts = counts, genes = rownames(counts)) # Again, move raw counts into x
x$samples$group <- group
keep.exprs <- rowSums(cpm > 1) >= 15
x <- x[keep.exprs,, keep.lib.sizes = FALSE]
dim(x)

lcpm <- cpm(x, log = T)

plot(
  density(lcpm[, 1]),
  col = col[1],
  ylim = c(0, 0.30),
  las = 2,
  main = "B. Filtered data",
  xlab = "Log-cpm"
)
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[, i])
  lines(den$x, den$y, col = col[i])
}

#================== 3. Normalize counts ================================================

# To properly normalize for library size we use TMM normalization, as discussed in the lectures.
x <- calcNormFactors(x, method = "TMM")
# If we look at the normalization factors, we see that the effect of TMM normalization was
# quite mild, which might not be too surprising with similar library sizes for each sample
x$samples

# Check after normalization:
par(mfrow = c(1, 2))
lcpm2 <- cpm(x, log = TRUE)
boxplot(lcpm2,  las = 2, col = col, main = "")
title(main = "A. normalised", ylab = "Log-cpm")


#================== 4. Unsupervised clustering ================================================

# Good to see if there are any outliers in the data:
par(mfrow = c(1, 2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels = group, col = col.group)
title(main = "A. Dimensions 1&2")
plotMDS(lcpm, labels = group, col = col.group, dim = c(3, 4))
title(main = "B. Dimensions 3&4")


# Look interactively with Glimma
glMDSPlot(lcpm, labels = group, groups = x$samples[, 2], folder = "scratch/glimma",
          launch = T)


#======================= 5. Pairwise DE analysis ==============================================
# We can analyze using pairwise contrasts, or using a generalized linear model. Here is
# demonstrated how to do pairwise.
# Using the grouping, we can now do an additional correction for the mean-variance correlation 
z <- estimateDisp(x)
# We can plot the mean-variance correlation
par(mfrow = c(1, 1))
plotBCV(z)
# The tagwise (individual genes) displacement from the common variance will be taken into
# account in remaining statistical framework.

# Regardless, we can now get lists of differentially expressed genes for each comparison, using
# the labels specified in group (let's have a look again what we specified:)
group
de2 <- exactTest(z, pair = c('ref', 'temp')) # Reference first!
# This will compare 'ref' with 'temp'. Show top 20 results:
topTags(de2, n = 20)

# Construct a table with all genes (n = Inf)
tt <- topTags(de2, n = Inf)

# Write CSV file if you prefer to inspect in e.g. Excel. Note that we're
# writing it to the results folder.
write.csv(tt, file = "results/ref_temp_DEgenes.csv", row.names = F)

# Make a smear plot, highlighting as DE those genes that have FDR < 0.01
plotSmear(x, pair = c('ref', 'pH'), de.tags = rownames(
  tt$table)[tt$table$FDR < 0.01 & abs(tt$table$logFC) > 1])


# Make histogram of p-values, add vertical line at p = 0.05
hist(tt$table$PValue, breaks = 100)
abline(v = 0.05, col = 'red')

# Make volcano plot
volcanoData <- cbind(tt$table$logFC, -log10(tt$table$FDR))
colnames(volcanoData) <- c("logFC", "-log10Pval")
plot(volcanoData)
abline(v = c(-2, 2), col = 'red')
abline(h = 3, col = 'red')

# Adjust the code above to make tables and CSV files for all comparisons that you are
# interested in.

# Also an option to make interactive graph with glimma:
de2 <- exactTest(z) # Reference first!
dt <- decideTests(de2)
glMDPlot(de2, status = dt, counts = x$counts, groups = x$samples$group, transform = TRUE,
         folder = "scratch/glimma")

#================== 6. DE analysis using linear model===========================================
# Specify design matrix as detailed in lecture
design <- model.matrix(~group)
colnames(design) <- gsub("group", '', colnames(design))
design

# Again calculate dispersion, now further informed by design matrix
z <- estimateDisp(z, design)

# Fit data to linear model
mfit <- glmQLFit(z, design)
# Inspect the coefficients that are found
head(mfit$coefficients)

# Perform DE analysis per factor
de <- glmLRT(mfit, coef = 2)
topTags(de)

# Or all de <- glmLRT(fit, coef=c(2,3,4,5))
de <- glmLRT(mfit, coef = c(2,3,4,5))
topTags(de)
tt <- topTags(de, n = Inf)
write.csv(tt, file = "results/glm_DEgenes.csv", row.names = F)

# Make clustering plot of DE genes only (similar as shown for proteomics data)
# First get a list of genes that are strongly differentially expressed upon the conditions
genes <- tt$table$genes[tt$table$FDR < 0.001]
plotMDS(lcpm[genes, ], labels = group, col = col.group)

#================== 7. Gene-set analysis =========================================
# Instead of using functions such as bioMart to import GO terms, you can also
# make your own list of gene-sets and load them as CSV file
#Note that gene name format has to be the same...for this file, I had to find and replace all _ with nothing to
#make GO list and counts in the same format.

mapGO <- read.csv('W29_CLIB122_GOTerm.csv')

#lots of work went into this, see 'Ortholog and GO Term K.marx.txt' and 'unpivot.txt'
#basically, DMKU gene names and GO terms were matched to CBS6556 gene names and rearranged for use in Piano

mapGO[1:10,]

# Load gene-sets for usage by Piano
myGsc <- loadGSC(mapGO)

head(mapGO)

# GSA is always on one statistic, so for instance one comparison between
# 2 conditions. Here we extract one of the coefficients of the GLM, and its
# associated P-value. As detailed in the lecture, the coefficients of the GLM
# represents the log2FC caused the the associated factor.

# Here we look at the second factor, which is temp
head(mfit$coefficients)

de <- glmLRT(mfit, coef = 2) # Fit linear model to extract the effect of temp
tt <- topTags(de, n = Inf) # Make a table with P-values and log2FC that will
# be used by Piano
#For pairwise, we called the DE dataset de2, so tt<- topTags(de2, n=Inf)

# Extract necessary data from table
Pval <- tt$table$PValue
names(Pval) <- tt$table$genes
FC <- tt$table$logFC
names(FC) <- tt$table$genes

# Run GSA with the extracted P-value, FC, the previously loaded gene-set
# collection myGsc, and we set a size limit of 10, 400.
gsaRes <- runGSA(Pval, FC, gsc = myGsc, gsSizeLim = c(3, 1000))

#For the purposes of comparison, I made 1 to 10000, 3 to 1000, and 10 to 400. GO numbers were 2188, 826, and 179 for K.marx
# For Y.lipo 1-10000 (), 3-1000 (702 gene sets)
# For Sac Cer 3-1000 (1390 gene sets)
#3+ seems the most appropriate

#To save the entire table
GSAsummaryTable(gsaRes, save = TRUE, file = 'low_pH_GSA3-1000.txt')

# We can show the results in a heatmap
GSAheatmap(gsaRes, cutoff = 5, cellnote = "rank", columnnames = "abbr", ncharLabel = "10")

# this takes a minute or two to generate
# Sometimes R comes with an error message regarding plotting,
# for instance "Error in plot.new() : figure margins too large". If the
# picture is plotted as intended, you can just ignore these messages.
# If the picture is not shown, it can help to run the following command
# (potentially multiple times, until it returns "null device 1":

# dev.off()


# NADH oxidation might be interesting, let's have a closer look.
# We can find out which genes are part of this geneset:
# (make sure you use the correct quotation mark around `GO term`, this
# is needed due to spaces in the GO terms.)
genes <- gsaRes$gsc$`arginine biosynthetic process`
# Extract the gene level data (note that this is for the anox coefficient)
tt[genes,]

# Maybe we want to see transcript levels of these genes for all samples.
# We'll use the normalized counts for this:
lcpm <- cpm(x)
boxData <- data.frame(lcpm[genes,], check.names = F)

# We make boxplots with ggplot2, which requires a specific data format (compare boxData
# before and after the following commands)
boxData <- cbind(genes = rownames(boxData), boxData) # Add genes as separate column
boxData <- gather(boxData, "sample", "count", 2:11) # Reorganize structure for ggplot2
map <- setNames(group, colnames(lcpm)) # Replicates should have same sample name. Map them here,...
boxData$sample <- map[boxData$sample] # and replace replicate name with sample name here.

# Make boxplots with ggplot2
ggplot(boxData, aes(sample, count)) +
  geom_boxplot(aes(colour = sample), size = 1) +
  facet_wrap( ~ genes, scales = "free_y") +  theme_bw() +
  labs(title = "arginine biosynthetic process", x = "", y = "Normalized counts")

# We can also connect the GO terms in a network plot
networkPlot(gsaRes, 'distinct', 'down',
            adjusted = T, significance = 0.01)

# Or show the results in a table
gsaTable <- GSAsummaryTable(gsaRes)
View(gsaTable)

# Run hypergeometric GSA, where we manually select a set of genes. Here, we
# select the most DE genes (FDR < 0.01) in high temperature (coef 5).
de <- glmLRT(mfit, coef = 5)
tt <- topTags(de, n = Inf)

hypGsaRes <- runGSAhyper(genes = tt$table$genes, pvalues = tt$table$FDR,
                         pcutoff = 0.01, gsc = myGsc)

networkPlot(hypGsaRes, 'non')

# Run consensus GSA. Let's speed up the process, we will only do 1000
# permutations when calculating the gene-set signficance (default is 10,000).
# To speed things up, we can check how many CPUs your computer has, and
# try to use them all
library(parallel)
cores <- as.numeric(detectCores())
# This is the nubmer of cores available:
cores
# Okay, back to consensus GSA:

gsaRes1 <- runGSA(Pval, directions = FC, geneSetStat = "mean", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes2 <- runGSA(Pval, directions = FC, geneSetStat = "median", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes3 <- runGSA(Pval, directions = FC, geneSetStat = "sum", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes4 <- runGSA(FC, geneSetStat = "maxmean", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes5 <- runGSA(Pval, directions = FC, geneSetStat = "fisher", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes6 <- runGSA(Pval, directions = FC, geneSetStat = "stouffer", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes7 <- runGSA(Pval, directions = FC, geneSetStat = "tailStrength", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes8 <- runGSA(FC, geneSetStat = "gsea", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)
gsaRes9 <- runGSA(FC, geneSetStat = "page", gsc = myGsc,
                  nPerm = 1000, gsSizeLim = c(10, 800), ncpus = cores)

# Once we have several GSA results, we can combine them and make a heatmap
resList <- list(gsaRes1, gsaRes2, gsaRes3, gsaRes4, gsaRes5, gsaRes6,
                gsaRes7, gsaRes8, gsaRes9)
names(resList) <- c("mean", "median", "sum", "maxmean", "fisher", "stouffer",
                    "tailStrength", "gsea", "page")
ch <- consensusHeatmap(resList, cutoff = 5, method = "mean")

# Heatmap creation in R
library (ade4)
library (gplots)
library(scatterplot3d)
library (made4)

# Aaron's method
heatplot(lcpm, margins=c(10,9), dend="row", distfun="euclidean", main="Y.lipo RNA-seq", runSampleTree="TRUE")