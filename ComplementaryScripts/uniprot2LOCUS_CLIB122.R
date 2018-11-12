#Open file to convert
setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Orthologs')
Yali_Ortho <- read.delim('Yali_annotation_uniprot.txt', header = TRUE, sep = "\t",stringsAsFactors=FALSE)
colnames(Yali_Ortho) <- c('OG_ID','CLIB122','W29')
#Open equivalences file
setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Databases/Uniprot')
YaliUni <- read.delim('uniprot-Yli_clib122.txt', header = TRUE, sep = "\t",stringsAsFactors=FALSE)
YaliUni <- YaliUni[,1:4]
YaliUni <- YaliUni[YaliUni[,4]>0,]

matches <- match(Yali_Ortho$CLIB122,YaliUni$Entry)
#Get indexes for converted elements
Yali_Ortho <- Yali_Ortho[!is.na(matches),]
converted <- matches[!is.na(matches)]
geneIDs   <- YaliUni[converted,4]

#Assign gene IDs for CLIB122 based on uniprot code
Yali_Ortho$CLIB122 <- YaliUni[converted,4]
setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Orthologs')
write.table(Yali_Ortho, "Yali_annotation.txt", sep="\t",row.names = FALSE)
