#multiOmics_Analysis_pipeline
#
#Function that performs differential expression analysis on proteomics absolute and relative) nd RNAseq datasets 
#for the datasets in the CHASSY project (fermentations for S. cerevisiae, K. marxianus and Y. lipolytica exposed
#reference, high temperature, low pH and osmotic stress conditions).
#
#The DE hits for all organisms and conditions are mapped to a list of 1:1:1 orthologous genes (from orthoFinder)
#to search for evolutionary conserved stress-adaptation responses at the transcript and protein levels.
#
#An intgrated table is also generated in which information of foldchanges at the transcript level together with 
#absolute proteomics levels [umol/g protein], Molecular weight of proteins, Sequence lenght and GO terms information
#is put together for each organism.
#
#This script will facilitate all analyses described above, the user only needs to 1) clone the repo and 
#2) change the directory name on Line 50 to reflect the location of your cloned directory
#
# Last modified: Ivan Domenzain. 2019-05-20
#

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
#Add the directory info where you cloned the OrthOmics Repository (the main directory), an example is shown
repoPath   <- '/Users/ivand/Documents/GitHub/OrthOmics'
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
source('mapDEgenesToOG.R')
mapDEgenesToOG(c('cpk','kma','yli'),pVal,log2FC,adjustedP,'RNA',repoPath)
