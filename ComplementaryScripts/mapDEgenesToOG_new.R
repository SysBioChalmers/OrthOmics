library(ggplot2)


# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
# Setting working directory to parent directory of ComplementaryScripts/
setwd("../")
# Assigning the working directory to the character object 'repoPath'
repoPath <- getwd()
dataPath <- paste(repoPath,'/RNA-seq/sce/Results/DE_log2FC_1_FDR_0.01',sep='')
# First, let's convert the uniprot IDs in the orthoGroups.txt fie (provided by Philip) to s2888c gene IDs
conditions <- c('HiT','LpH','Osm')
#Read orthoGroups file
#For j-th organism get the data for the genes that were DE exclusively in the i-th condition
filename    <- paste(repoPath,'/Databases/orthoGroups.txt',sep='')
OG_df       <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
filename    <- paste(repoPath,'/Databases/Uniprot/uniprot_sce.txt',sep='')
uniprot     <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
#map genes id OG df to uniprot df
indexesOG   <- match(OG_df$gene,uniprot$Entry)
#Translate IDs
OrthoGroups <- OG_df
OrthoGroups$Stratum[is.na(OrthoGroups$Stratum)] <- 'NaN'
OrthoGroups$Stratum[OrthoGroups$Stratum==5] <- 4
OrthoGroups$gene <- uniprot$Gene.names...ordered.locus..[indexesOG]
#Get information about the Sce genes grouping
gene_groups <- unique(OrthoGroups$Stratum)
gene_groups <- gene_groups[order(gene_groups)]
group_totals <- c()
for (i in 1:length(gene_groups)){
  group_totals <- c(group_totals,sum(is.element(OrthoGroups$Stratum,gene_groups[i])))
}
#Now let's go through all conditions and map DE genes to OrthoGroups
resultsDF <- group_totals
for (i in 1:length(conditions)){
  filename <- paste(dataPath,'/sce_RNA_ref_',conditions[i],'.csv',sep='')
  DEdata   <- read.csv(filename,stringsAsFactors = FALSE)
  allGenes <- DEdata$genes
  DE_genes <- DEdata$genes[abs(DEdata$logFC)>=1 & DEdata$FDR<=0.01]
  DE_ratio <- length(DE_genes)/length(DEdata$genes)
  #Map all genes in dataset to OG grouping list
  group_totals  <- c()
  group_DE      <- c()
  group_DE_perc <- c()
  colors <- c('blue','red','grey','yellow','purple','cyan','brown')
  alphaVal <- 0.5
  for (j in 1:length(gene_groups)){
    cond_genes_groups <- OrthoGroups$gene[which(OrthoGroups$Stratum==gene_groups[j])]
    group_totals  <- c(group_totals,length(cond_genes_groups))
    group_DE      <- c(group_DE,sum(is.element(DE_genes,cond_genes_groups)))
    group_DE_perc <- c(group_DE_perc,(group_DE[j]/group_totals[j])/DE_ratio)
  }
  dt<-data.frame(gene_groups,group_totals,group_DE,group_DE_perc)
  dt$categories <- c('0','1','2','3','4','NA')
  fontSize <- 18
  p <- ggplot(data=dt, aes(x=categories,y=(group_DE_perc))) + geom_bar(stat='identity')+ 
     theme_bw(base_size = 2*fontSize)+xlab('Gene groups') +
    ylab('fold_DE_enrichment')#+ scale_y_continuous(breaks = (seq(0,5,by = 0.5)))
  plotTitle <- paste(repoPath,'/newMapping/sce_',conditions[i],'_AllDE.png',sep='')
  png(plotTitle,width = 600, height = 600)
  plot(p)
  dev.off()
  #+ xlab(xLabel) +
}