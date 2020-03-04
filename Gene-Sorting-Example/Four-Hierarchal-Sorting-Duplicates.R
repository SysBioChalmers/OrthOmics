#Load Duplicate Gene Sorting Results and All protein coding genes sheet

setwd('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting')
GroupI_Duplicates  <- read.csv("GroupIDUPEs_SCE.csv", header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
GroupII_Duplicates  <- read.csv("GroupIIDUPEs_SCE.csv", header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
GroupIII_Duplicates  <- read.csv("GroupIIIDUPEs_SCE.csv", header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
GroupIV_Duplicates  <- read.csv("GroupIVDUPEs_SCE.csv", header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
gene_2_prot <- read.csv("s288c_gene_prot.csv", header = TRUE)
AllDuplicates <- read.csv('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/GroupV/Duplicates_to_Exclude.csv', header = TRUE,  sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
AllDuplicates <- AllDuplicates[c(1,2)]
TransposonsToExclude <- read.csv("transposons.csv", header = TRUE)
WholeGenomeDupes <- read.csv('WGD-Byrne2005.csv', header = TRUE,  sep = ",",stringsAsFactors=FALSE, na.strings = "#N/A")
Three_Plus_CopyGenes <- read.csv('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/GroupV/3+Copy_Genes.csv', header = TRUE,  sep = ",",stringsAsFactors=FALSE, na.strings = "#N/A")
Three_Plus_CopyGenes <- Three_Plus_CopyGenes[c(1,2)]

#First we need to remove some WGD genes that have reduplicated (exist as 3+ copy genes)
WholeGenomeDupes <- (WholeGenomeDupes %>% anti_join(Three_Plus_CopyGenes))

#Duplicate sorting is bottom-UP - starts with groupV and tests which duplicates can be traced backwards
AllDuplicates <- (AllDuplicates %>% anti_join(TransposonsToExclude))
AllDuplicates <- (AllDuplicates %>% anti_join(WholeGenomeDupes))

GroupV <- (AllDuplicates %>% anti_join(GroupIV_Duplicates))
namevector <- c("GroupV")
GroupV[, namevector] <- "GroupV"
names(GroupV)[names(GroupV) == "GroupV"] <- "Group"

GroupIV <- (AllDuplicates %>% anti_join(GroupV))
GroupIV <- (GroupIV %>% anti_join(GroupIII_Duplicates))
namevector <- c("GroupIV")
GroupIV[, namevector] <- "GroupIV"
names(GroupIV)[names(GroupIV) == "GroupIV"] <- "Group"

GroupIII <- (AllDuplicates %>% anti_join(GroupV))
GroupIII <- (GroupIII %>% anti_join(GroupIV))
GroupIII <- (GroupIII %>% anti_join(GroupII_Duplicates))
namevector <- c("GroupIII")
GroupIII[, namevector] <- "GroupIII"
names(GroupIII)[names(GroupIII) == "GroupIII"] <- "Group"

GroupII <- (AllDuplicates %>% anti_join(GroupV))
GroupII <- (GroupII %>% anti_join(GroupIV))
GroupII <- (GroupII %>% anti_join(GroupIII))
GroupII <- (GroupII %>% anti_join(GroupI_Duplicates))
namevector <- c("GroupII")
GroupII[, namevector] <- "GroupII"
names(GroupII)[names(GroupII) == "GroupII"] <- "Group"

GroupI <- (AllDuplicates %>% anti_join(GroupV))
GroupI <- (GroupI %>% anti_join(GroupIV))
GroupI <- (GroupI %>% anti_join(GroupIII))
GroupI <- (GroupI %>% anti_join(GroupII))
namevector <- c("GroupI")
GroupI[, namevector] <- "GroupI"
names(GroupI)[names(GroupI) == "GroupI"] <- "Group"

#Assemble hierarchy for non-WGD Duplicates
SCESortedDuplicates <- rbind(GroupI, GroupII, GroupIII, GroupIV, GroupV)

#Add protein-coding WGD genes to the list
WholeGenomeDupes <- WholeGenomeDupes[complete.cases(WholeGenomeDupes[ ,2 ]), ]
namevector <- c("WGD")
WholeGenomeDupes[, namevector] <- "WGD"
names(WholeGenomeDupes)[names(WholeGenomeDupes) == "WGD"] <- "Group"
SCESortedDuplicates <- rbind(GroupI, GroupII, GroupIII, WholeGenomeDupes, GroupIV, GroupV)

write.csv(SCESortedDuplicates, file="Sorted_MultiCopy_Genes_SCE.csv", row.names = FALSE)

SortedSingles <- read.csv("Sorted_SingleCopy_Genes_SCE.csv", header = TRUE)

SortedAll <- rbind(SCESortedDuplicates, SortedSingles)
write.csv(SortedAll, file="All_Genes_Sorted_SCE.csv", row.names = FALSE)
