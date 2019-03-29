plotVennDiagram <- function(elementsList,categories,colorValues,intLabSize,ellipses,scaleF){
library(VennDiagram)
  nargin <- length(as.list(match.call())) -1
  if (nargin < 6){scaleF = FALSE}
  euler <- TRUE
  overlap <- list()  
    if (ellipses  >2){
    area   <- c()
    shared <- c()
    for (i in 1:3){
      if (i==3){j<-1}
      else {j<-i+1}
      area[i]      <- length(elementsList[[i]])
      overlap[[i]] <- intersect(elementsList[[i]],elementsList[[j]])
      shared[i]    <- length(overlap[[i]])
    }
    #Intersect 1,2,3
    overlap[[4]] <- intersect(intersect(elementsList[[1]],elementsList[[2]]),elementsList[[3]])
    shared[4] <- length(overlap[[4]])
    
    if (ellipses <4){
    venn.plot<-draw.triple.venn(area[1], area[2],area[3], 
                                shared[1],shared[2],shared[3],shared[4],
                                category = categories,scaled = scaleF,euler.d = euler,fill = colorValues, 
                                lty = c(rep("solid",3)), cex = intLabSize, cat.cex = 2.5)}
    else{
      area[4] <- length(elementsList[[4]])
      #Intersect 2,4
      shared[5] <- length(intersect(elementsList[[2]],elementsList[[4]]))
      #Intersect 3,4
      shared[6] <- length(intersect(elementsList[[3]],elementsList[[4]]))
      #Intersect 4,1
      shared[7] <- length(intersect(elementsList[[1]],elementsList[[4]]))
      #Intersect 1,2,4
      shared[8] <- length(intersect(intersect(elementsList[[1]],elementsList[[2]]),elementsList[[4]]))
      #Intersect 1,3,4
      shared[9] <- length(intersect(intersect(elementsList[[1]],elementsList[[3]]),elementsList[[4]]))
      #Intersect 2,3,4
      shared[10] <- length(intersect(intersect(elementsList[[2]],elementsList[[3]]),elementsList[[4]]))
      #Intersect 1,2,3,4
      shared[11] <- length(intersect(intersect(elementsList[[1]],elementsList[[2]]),intersect(elementsList[[3]],elementsList[[4]]))) 

      venn.plot <- draw.quad.venn(area[1], area[2],area[3],area[4], 
                                  shared[1],shared[3],shared[7],shared[2],shared[5],
                                  shared[6],shared[4],shared[8],shared[9],shared[10],shared[11],
                                  #n12, n13,n14, n23, n24,
                                  #n34, n123, n124, n134, n234, n1234, 
                                  category = categories,scaled = scaleF,fill = colorValues, 
                                  lty = c(rep("solid",4)), cex = intLabSize, cat.cex = 2.5,euler.d = euler)
      
    }
  }
  if (ellipses  == 2){
    area   <- c()
    shared <- c()
    for (i in 1:2){
      area[i]      <- length(elementsList[[i]])
      overlap[[i]] <- elementsList[[i]]
    } 
    overlap[[3]] <- intersect(elementsList[[1]],elementsList[[2]])
    #Intersect 1,2
    shared  <- length(overlap[[3]])
    venn.plot<-draw.pairwise.venn(area[1], area[2],shared,
                                category = categories,scaled = scaleF,fill = colorValues, 
                                lty = c(rep("solid",2)), cex = intLabSize, cat.cex = 2.5,euler.d = euler)
    

  }
  return(overlap)
}