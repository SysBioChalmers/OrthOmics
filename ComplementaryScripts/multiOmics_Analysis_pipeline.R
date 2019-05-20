install.packages('VennDiagram')
install.packages("rlist")
library("rlist")
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
organisms    <- c('sce','kma','yli')
#Indicate the dataset that should be used foer DE analysis
dataSource <- 'XIC' #1 for XIC or NSAF, 2 for Scounts or iBAQ and 3 for merged datasets
#Define DE thresholds
pVal       <- 0.01
logPval    <- abs(log10(pVal))
log2FC     <- 1
adjustedP  <- TRUE
#Should the initial dataset be normalized by MW of proteins
normByMW   <- TRUE
#Filter type for determination of present and non-noisy proteins in the dataset (TRUE if filter
#criteria should be applied to all conditions, FALSE if just the reference is desired to be 
#filtered)
stringent  <- TRUE
#Normalization method for DE analysis
normMethod <- 'TMM'
repoPath   <- '/Users/ivand/Documents/GitHub/CHASSY_multiOmics_Analysis'
#Internal functions path
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')

#================== Data analysis pipeline ====================================
for (organism in organisms){
  cat(paste("Analyzing data for: ", organism,'\n',sep=""))
  #================== 1. Analyze Transcriptomics Data ====================================
  setwd(scriptsPath)
  source('RNAseqAnalysis.R')
  if (all(organism=='sce')){orgID <- 'cpk'} else{orgID <- organism}
  cat(paste("Analyzing RNAseq data for: ", organism,'\n',sep=""))
  RNAseqAnalysis(orgID,stringent,normMethod,0,logPval,log2FC,adjustedP,repoPath)
  #================== 2. Analyze proteomics Data ====================================
  setwd(scriptsPath)
  source('proteomics_Analysis.R')
  cat(paste("Analyzing proteomics data for: ", organism,'\n',sep=""))
  proteomics_Analysis(organism,dataSource,normByMW,stringent,normMethod,logPval,log2FC,adjustedP,repoPath)
  #================== 3. Create integratative Omics table ==================
  setwd(scriptsPath)
  source('createIntegratedTable.R')
  createIntegratedTable(orgID,pVal,log2FC,adjustedP,FALSE)
  cat("\014") 
}
#================== 3. Map DE genes to 1:1:1 orthologous genes list ==================
setwd(scriptsPath)
source('mapDEgenesToOG.R')cat(paste("Analyzing proteomics data for: ", organism,'\n',sep=""))
mapDEgenesToOG(c('cpk','kma','yli'),pVal,log2FC,adjustedP,'RNA',repoPath)
