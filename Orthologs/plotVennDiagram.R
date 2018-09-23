plotVennDiagram <- function(elementsList,categories,colorValues,intLabSize,ellipses){
library(VennDiagram)
    if (ellipses  == 3){
    area   <- c()
    shared <- c()
    for (i in 1:3){
      if (i==3){j<-1}
      else {j<-i+1}
      area[i]   <- length(elementsList[[i]])
      shared[i] <- length(intersect(elementsList[[i]],elementsList[[j]]))
    }
    #Intersect 1,2,3
    overlap   <-intersect(intersect(elementsList[[1]],elementsList[[2]]),elementsList[[3]])
    shared[4] <- length(overlap)
    overlap
    
    venn.plot<-draw.triple.venn(area[1], area[2],area[3], 
                                shared[1],shared[2],shared[3],shared[4],
                                category = categories,scaled = F,fill = colorValues, 
                                lty = c(rep("solid",3)), cex = intLabSize, cat.cex = 2)
  }
  if (ellipses  == 2){
    area   <- c()
    shared <- c()
    for (i in 1:2){
      area[i]   <- length(elementsList[[i]])
    } 
    #Intersect 1,2,3
    overlap <-intersect(elementsList[[1]],elementsList[[2]])
    shared  <- length(overlap)
    venn.plot<-draw.pairwise.venn(area[1], area[2],shared,
                                category = categories,scaled = F,fill = colorValues, 
                                lty = c(rep("solid",2)), cex = c(1,2,1), cat.cex = 2)
  }
  return(overlap)
}