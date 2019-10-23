#This script allows self-self Orthology searches to be quickly stripped down to a list of
#genes and proteins that are predicted to be duplicates in an organism.

#This results in two outputs that are needed in future steps
#1) A list of genes to exclude from singleCopy Analysis called Duplicates_to_Exclude.csv
#2) A list of duplicated genes in the same ortholog family which we will later sort to determine the relative
#timing of the duplication event (Duplications_Groups_For_Sorting.csv)

setwd('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/GroupV')
filename      <- 'Orthogroups.GeneCount.tsv'
DuplicateOGs  <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
filename      <- 'Orthogroups.tsv'
orthoGroups   <- read.delim(filename, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
endCol        <- ncol(orthoGroups)
newOGlist     <- c()
group_number  <- 'Duplicates_to_Exclude.csv' #Add the group number for the analysis
filename   <- '/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/s288c_gene_prot.csv'
gene_2_prot <- read.csv(filename, header = TRUE)


DuplicateOGs <- DuplicateOGs[apply(DuplicateOGs[c(2)],1,function(z) any(z!=0)),] #Organism of interest in column 2 - Remove orthogroups with 0 values for the organism of interest (these are orthologs between the other organisms but not the query)
DuplicateOGs <- DuplicateOGs[apply(DuplicateOGs[c(2)],1,function(z) any(z>=2)),] # grab duplicates from self-self orthology search



#Map each  OG to orthoGroups data frame
for (i in 1:nrow(DuplicateOGs)){
  OGid    <- DuplicateOGs[i,1] 
  index   <- which(orthoGroups[,1]==OGid)
  rowData <- orthoGroups[index,2:endCol]
  newRow  <- OGid
  
  
#Split  protein ID's for each organism
  for (j in 1:(endCol-1)){
    cellText  <- rowData[[j]]

    newRow    <- cbind(newRow,cellText)
  }
  newOGlist <- rbind(newOGlist,newRow)
}
colnames(newOGlist) <- colnames(orthoGroups)
return(newOGlist)



newOGlistP <- newOGlist [,c(1:2)]
newOGlistP2 <- cSplit(newOGlistP, 'A_s288c', sep=",", type.convert=FALSE)

library(stringr)
#Keep only the uniprot protein ID

test3 <- lapply(newOGlistP2, word, 1)

test3 <- as.data.frame(test3)
rownames(test3) <- test3$Orthogroup
test3 <- test3[,-1]
#test3t <- t(test3)
#test3t <- colnames(test3t, test3t$Orthogroup)
#test3t <- as.data.frame(test3t, col.names = test3t$Orthogroup)
#test3t['V1']
#colnames(test3t) <- test3t$Orthogroup

#test3 <- test3[,-1]
#test3 <- toString(test3, width=NULL)
write.table(test3, file="Duplications_Groups_For_Sorting.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ", ")

#XXX <- read.csv ("Duplications_Groups_For_Sorting.csv", sep = " ")
#test4 <- lapply(test3, function(x) x[!is.na(x)])
#test4 <- as.matrix(test4)
#write(test4, file = "Duplications_As_Groups.txt")

#Pull all uniprot IDs into a linear list                    
newOGlistQ <- newOGlist [,c(2)]
newOGsimple <- toString(newOGlistQ, width=NULL)
test <- as.list(strsplit(newOGsimple, ","))
num.el <- sapply(test, length)
newOGsimple <- cbind(unlist(test), rep(1:length(test), num.el))
colnames(newOGsimple) <- c("protein","sortingGroup")
newOGsimple <- as.data.frame(newOGsimple)
library(splitstackshape)
newOGsimple2 <- cSplit(newOGsimple, 'protein', sep=" ", type.convert=FALSE)
newOGsimple2$sortingGroup <- gsub('1', 'Self-Self-Duplicate', newOGsimple$sortingGroup)
newOGsimple2 <- newOGsimple2[,c(1,2)]
newOGsimple2$protein <- gsub('\\s+', '', newOGsimple2$protein_01)


m <- merge(newOGsimple2, gene_2_prot, by.x = "protein", by.y = "protein", all.x= TRUE)
m  <- m [,c(1,4,2)]

#n <- merge(m, test3, by.x = c("protein"), by.y = c("A_s288c_01"), all.x= FALSE, all.y=FALSE)

write.csv(m, file = group_number, row.names = FALSE)


