First_ID_selfSelf_duplicates <- function(repoPath){
#This script allows self-self Orthology searches to be quickly stripped down to a list of
#genes and proteins that are predicted to be duplicates in an organism.

#This results in two outputs that are needed in future steps

#1) A list of duplicated genes in the same ortholog family which we will later sort to determine the relative
#timing of the duplication event (Duplications_Groups_For_Sorting.csv)
#2) A list of genes that exist as 3 or more copies for other analyses
#3) A list of genes to exclude from singleCopy Analysis called Duplicates_to_Exclude.csv

library(splitstackshape)
library(stringr)
repoPath <- '/Users/ivand/Documents/GitHub/OrthOmics'
#Load Necessary Data from the Self-Self orthology search
groupingPath <- paste(repoPath,'/Gene-Sorting-Example/S.cerevisiae_GeneSorting',sep='')
setwd(groupingPath)
filename     <- paste(groupingPath,'/GroupV/Orthogroups.GeneCount.tsv',sep='')
DuplicateOGs <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
filename     <- paste(groupingPath,'/GroupV/Orthogroups.tsv',sep='')
orthoGroups  <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
endCol       <- ncol(orthoGroups)
group_number <- 'Duplicates_to_Exclude.csv' #Add the group number for the analysis
filename     <- paste(groupingPath,'/s288c_gene_prot.csv',sep='')
gene_2_prot  <- read.csv(filename, header = TRUE)
newOGlist    <- c()
setwd(paste(groupingPath,'/GroupV',sep=''))
#Keep only self-self duplicates (i.e. 2+ copy genes)
rownames(DuplicateOGs) <- DuplicateOGs[,1]
DuplicateOGs           <- DuplicateOGs[,2:endCol]
#Organism of interest in column 1 - Remove orthogroups with 0 values for the organism 
#of interest (these are orthologs between the other organisms but not the query)
DuplicateOGs <- DuplicateOGs[apply(DuplicateOGs[c(1)],1,function(z) any(z!=0)),] 
# grab duplicates from self-self orthology search
DuplicateOGs <- DuplicateOGs[apply(DuplicateOGs[c(1)],1,function(z) any(z>=2)),] 
#Map each  Gene/protein to orthoGroups data frame
matches   <- match(rownames(DuplicateOGs),orthoGroups[,1])
newOGlist <- orthoGroups[matches,]
colnames(newOGlist) <- colnames(orthoGroups)
#return(newOGlist)
newOGlistP <- newOGlist [,c(1:2)]
newOGlistP2 <- cSplit(newOGlistP, 'A_s288c', sep=",", type.convert=FALSE)
#1) Create a list of duplicates for sorting later on in the format that will match orthofinder outputs
#Keep the protein IDs inside of their respective orthogroups for future sorting
MultiCopyGroups <- lapply(newOGlistP2, word, 1)
MultiCopyGroups <- as.data.frame(MultiCopyGroups)
rownames(MultiCopyGroups) <- MultiCopyGroups$Orthogroup
MultiCopyGroups <- MultiCopyGroups[,-1]
write.table(MultiCopyGroups, file="Duplications_Groups_For_Sorting.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ", ", na="")
#2) Create a list of genes that have 3+copies in the organism of interest - This helps for downstream analyses
Three_Plus_Copy_Genes <- MultiCopyGroups[complete.cases(MultiCopyGroups[ ,3 ]), ]
Three_Plus_Copy_Genes <- data.frame(matrix(unlist(Three_Plus_Copy_Genes), byrow=T), stringsAsFactors=FALSE)
namevector <- c("3+CopyGenes")
Three_Plus_Copy_Genes[, namevector] <- "3+CopyGenes"
names(Three_Plus_Copy_Genes)[names(Three_Plus_Copy_Genes) == "3+CopyGenes"] <- "Group"
names(Three_Plus_Copy_Genes) [1] <- paste("protein")
Three_Plus_Copy_Genes <- Three_Plus_Copy_Genes[complete.cases(Three_Plus_Copy_Genes[ ,1 ]), ]
Three_Plus_Copy_Genes <- merge(Three_Plus_Copy_Genes, gene_2_prot, by.x = "protein", by.y = "protein", all.x= TRUE)
Three_Plus_Copy_Genes <- Three_Plus_Copy_Genes[,c(1,3,2)]
write.table(Three_Plus_Copy_Genes, file="3+Copy_Genes.csv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
#3) Create a list of duplicated genes to exclude from singlegene sorting
#Pull all uniprot IDs into a linear list                    
newOGlistQ <- newOGlist [,c(2)]
newOGsimple <- toString(newOGlistQ, width=NULL)
test <- as.list(strsplit(newOGsimple, ","))
num.el <- sapply(test, length)
newOGsimple <- cbind(unlist(test), rep(1:length(test), num.el))
colnames(newOGsimple) <- c("protein","sortingGroup")
newOGsimple <- as.data.frame(newOGsimple)
newOGsimple2 <- cSplit(newOGsimple, 'protein', sep=" ", type.convert=FALSE)
newOGsimple2$sortingGroup <- gsub('1', 'Self-Self-Duplicate', newOGsimple$sortingGroup)
newOGsimple2 <- newOGsimple2[,c(1,2)]
newOGsimple2$protein <- gsub('\\s+', '', newOGsimple2$protein_01)


m <- merge(newOGsimple2, gene_2_prot, by.x = "protein", by.y = "protein", all.x= TRUE)
m  <- m [,c(1,4,2)]

write.csv(m, file = group_number, row.names = FALSE)
}

