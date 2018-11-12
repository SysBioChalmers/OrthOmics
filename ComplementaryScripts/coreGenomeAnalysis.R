#Load a FC data file
setwd('/Users/ivand/Documents/GitHub/CHASSY-Multi-Omics-Analyisis/Orthologs')
source('plotVennDiagram.r')

matrix           <- Most_changed_genes[,2:9]
rownames(matrix) <- Most_changed_genes$`OG Name`
nConditions      <- c(3,2,3)
conditions       <- c('HiT','LpH','Osm')
nOrganisms       <- 3
organisms    <- c('Sce','Yli','Kma')
colorValues  <- colorValues <- c("blue", "yellow", "red")
#Get a transformed matrix with the sign of the elements
sign.matrix   <- sign(matrix)

#Get the genes that were up/down-regulated in all conditions
AllDown       <- rownames(matrix)[rowSums(sign.matrix < 0) == ncol(sign.matrix)]
AllUp         <- rownames(matrix)[rowSums(sign.matrix > 0) == ncol(sign.matrix)]
#Now those that were up/down-regulated in at least 75% of the conditions  
mostDown      <- rownames(matrix)[rowSums(sign.matrix < 0) >= 0.75*ncol(sign.matrix)]
mostUp        <- rownames(matrix)[rowSums(sign.matrix > 0) >= 0.75*ncol(sign.matrix)]

#Now get a matrix That shows those genes up or down regulated for all the conditions for each of the organisms
orgDOWN <- c()
orgDOWN[[1]] <-  rownames(matrix)[rowSums(sign.matrix[,1:3] < 0) == nConditions[1]]
orgDOWN[[2]] <-  rownames(matrix)[rowSums(sign.matrix[,4:5] < 0) == nConditions[2]]
orgDOWN[[3]] <-  rownames(matrix)[rowSums(sign.matrix[,6:8] < 0) == nConditions[3]]
orgUP <- c()
orgUP[[1]] <-  rownames(matrix)[rowSums(sign.matrix[,1:3] > 0) == nConditions[1]]
orgUP[[2]] <-  rownames(matrix)[rowSums(sign.matrix[,4:5] > 0) == nConditions[2]]
orgUP[[3]] <-  rownames(matrix)[rowSums(sign.matrix[,6:8] > 0) == nConditions[3]]

#Get a venn diagram with the genes that were DOWN regulated for all conditions
intLabSize   <- c(rep(1,7))
intLabSize[5]<- 2
plotVennDiagram(orgDOWN,organisms,colorValues,intLabSize,3)
#Get a venn diagram with the genes that were UP regulated for all conditions
intLabSize   <- c(rep(1,7))
intLabSize[2]<- 2
intLabSize[4]<- 2
intLabSize[6]<- 2
plotVennDiagram(orgUP,organisms,colorValues,intLabSize,3)

#Now get venn diagrams for up and down regulated genes per condition
overlapUp <- c()
overlapDown <- c()
for (i in 1:length(conditions)){

  if (i==1){condCols<- c(3,5,6)}
  if (i==2){condCols<- c(2,4,8)}
  if (i==3){
    condCols<- c(1,7)
    organisms <- c('Sce','Kma')
    colorValues  <- colorValues <- c("blue", "red")
  }
  ellipses <- length(condCols)
  #Extract columns for the condition of interest
  condMat     <- matrix[,condCols]
  sign.matrix <- sign(condMat)
  Down <- c()
  Up <- c()
  for (j in 1:length(condCols)){
    orgMat    <- sign.matrix[,j]
    Down[[j]]     <- rownames(condMat)[orgMat < 0]
    Up[[j]]       <- rownames(condMat)[orgMat > 0]
  }
  indexes <- plotVennDiagram(Down,organisms,colorValues,intLabSize,ellipses)
  overlapDown[[i]] <- rownames(matrix)[as.numeric(indexes)]
  indexes <- plotVennDiagram(Up,organisms,colorValues,intLabSize,ellipses)
  overlapUp[[i]] <- rownames(matrix)[as.numeric(indexes)]
}
#Get the evolutionary conserved common stress responses (HiT and LpH)
conservedResp_D <- intersect(overlapDown[[1]],overlapDown[[2]])
conservedResp_U <- intersect(overlapUp[[1]],overlapUp[[2]])