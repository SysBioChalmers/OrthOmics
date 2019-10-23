#Load Single Copy Gene Sorting Results and All protein coding genes sheet

setwd('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting')
GroupI_Singles  <- read.csv("GroupI-Singles.csv", header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
GroupII_Singles  <- read.csv("GroupII-Singles.csv", header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
GroupIII_Singles  <- read.csv("GroupIII-Singles.csv", header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
GroupIV_Singles  <- read.csv("GroupIV-Singles.csv", header = TRUE, sep = ",",stringsAsFactors=FALSE, na.strings = "NA")
gene_2_prot <- read.csv("s288c_gene_prot.csv", header = TRUE)
DuplicatesToExclude <- read.csv('/Users/doughty/Documents/GitHub/CHASSY_multiOmics_Analysis/Gene-Sorting-Example/S.cerevisiae_GeneSorting/GroupV/Duplicates_to_Exclude.csv', header = TRUE,  sep = ",",stringsAsFactors=FALSE, na.strings = "NA")

SingleCopyGenes <- (gene_2_prot %>% anti_join(DuplicatesToExclude))

#Begin hierarchal sorting
GroupI <- subset(SingleCopyGenes, SingleCopyGenes$protein %in% GroupI_Singles$protein)
namevector <- c("GroupI")
GroupI[, namevector] <- "GroupI"
names(GroupI)[names(GroupI) == "GroupI"] <- "Group"
subtractI <- (SingleCopyGenes %>% anti_join(GroupI))


GroupII <- subset(subtractI, subtractI$protein %in% GroupII_Singles$protein)
namevector <- c("GroupII")
GroupII[, namevector] <- "GroupII"
names(GroupII)[names(GroupII) == "GroupII"] <- "Group"
subtractII <- (subtractI %>% anti_join(GroupII))

GroupIII <- subset(subtractII, subtractII$protein %in% GroupIII_Singles$protein)
namevector <- c("GroupIII")
GroupIII[, namevector] <- "GroupIII"
names(GroupIII)[names(GroupIII) == "GroupIII"] <- "Group"
subtractIII <- subset(subtractII %>% anti_join(GroupIII))

GroupIV <- subset(subtractIII, subtractIII$protein %in% GroupIV_Singles$protein)
namevector <- c("GroupIV")
GroupIV[, namevector] <- "GroupIV"
names(GroupIV)[names(GroupIV) == "GroupIV"] <- "Group"
subtractIV <- subset(subtractIII %>% anti_join(GroupIV))

GroupV <- subtractIV
namevector <- c("GroupV")
GroupV[, namevector] <- "GroupV"
names(GroupV)[names(GroupV) == "GroupV"] <- "Group"

#Assemble hierarchy
SCESortedSingles <- rbind(GroupI, GroupII, GroupIII, GroupIV, GroupV)

write.csv(SCESortedSingles, file="Sorted_SingleCopy_Genes_SCE.csv", row.names = FALSE)
