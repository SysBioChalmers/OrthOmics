#AnalysisPipeline_mainScript
#
#Function that performs differential expression analysis on proteomics absolute and relative) nd RNAseq datasets 
#for the datasets in the CHASSY project (fermentations for S. cerevisiae, K. marxianus and Y. lipolytica exposed
#reference, high temperature, low pH and osmotic stress conditions).
#
#The DE hits for all organisms and conditions are mapped to a list of 1:1:1 orthologous genes (from orthoFinder)
#to search for evolutionary conserved stress-adaptation responses at the transcript and protein levels.
#
#An integrated table is also generated in which information of foldchanges at the transcript level together with 
#absolute proteomics levels [umol/g protein], Molecular weight of proteins, Sequence lenght and GO terms information
#is put together for each organism.
#
#This script will facilitate all analyses described above, the user only needs to 1) clone the repo
#
# Last modified: Ronan Harrington. 2020-04-27
#

# Changing working directory ===================================================

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

# Setting working directory to parent directory of ComplementaryScripts/
setwd("../")

# Installing required R packages ===============================================
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")}
if (!requireNamespace("VennDiagram", quietly = TRUE)){
  install.packages("VennDiagram")}
if (!requireNamespace("rlist", quietly = TRUE)){
  install.packages("rlist")}
if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")}
if (!requireNamespace("RColorBrewer", quietly = TRUE)){
  install.packages("RColorBrewer")}
if (!requireNamespace("ggbiplot", quietly = TRUE)){
  install.packages("ggbiplot")}
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
if (!requireNamespace("limma", quietly = TRUE)){
  BiocManager::install("limma")}
if (!requireNamespace("edgeR", quietly = TRUE)){
  BiocManager::install("edgeR")}

# Loading libraries ============================================================
library("rlist")
library(VennDiagram)
library(devtools)
library(ggplot2)
library(limma)
library(edgeR)
library(tidyverse) # Collection of useful R-packages
library(RColorBrewer)
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
# Assigning the working directory to the character object 'repoPath'
repoPath   <- getwd()
#Internal functions path
scriptsPath <- paste(repoPath,'/ComplementaryScripts',sep='')

#================== Data analysis pipeline ====================================
for (organism in organisms){
  cat(paste("Analyzing data for: ", organism,'\n',sep=""))
  #================== 1. Analyze Transcriptomics Data ====================================
  setwd(scriptsPath)
  source('RNAseqAnalysis.R')
  cat(paste("Analyzing RNAseq data for: ", organism,'\n',sep=""))
  RNAseqAnalysis(organism,normMethod,0,logPval,log2FC,adjustedP,repoPath)
  #================== 2. Analyze proteomics Data ====================================
  setwd(scriptsPath)
  source('proteomics_Analysis.R')
  cat(paste("Analyzing proteomics data for: ", organism,'\n',sep=""))
  proteomics_Analysis(organism,dataSource,normByMW,stringent,normMethod,logPval,log2FC,adjustedP,repoPath)
  #================== 3. Create integratative Omics table ==================
  setwd(scriptsPath)
  source('createIntegratedTable.R')
  createIntegratedTable(organism,pVal,log2FC,adjustedP,FALSE,repoPath)
  cat("\014") 
}
#================== 3. Map DE genes to 1:1:1 orthologous genes list ==================
setwd(scriptsPath)
source('mapDEgenesToOG.R')
mapDEgenesToOG(c('sce','kma','yli'),pVal,log2FC,adjustedP,'RNA',repoPath)
