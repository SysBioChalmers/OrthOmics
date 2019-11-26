plotVennDiagram <- function(elementsList,categories,colorValues,intLabSize,ellipses,scaleF){
library(VennDiagram)
  nargin <- length(as.list(match.call())) -1
  if (nargin < 6){scaleF = FALSE}
  euler <- TRUE
  shared12  <- length(intersect(elementsList[[1]],elementsList[[2]]))
  overlap12 <- intersect(elementsList[[1]],elementsList[[2]])
  if (ellipses  >2){
    area   <- c()
    shared <- c()
    for (i in 1:3){
      if (i==3){j<-1}
      else {j<-i+1}
      area[i]      <- length(elementsList[[i]])
    }
    shared23 <- length(intersect(elementsList[[2]],elementsList[[3]]))
    shared13 <- length(intersect(elementsList[[1]],elementsList[[3]]))
    overlap123 <- intersect(intersect(elementsList[[1]],elementsList[[2]]),elementsList[[3]])
    shared123 <- length(overlap123)
    if (ellipses >=4){
        area[4] <- length(elementsList[[4]])
        shared24 <- length(intersect(elementsList[[2]],elementsList[[4]]))
        shared34 <- length(intersect(elementsList[[3]],elementsList[[4]]))
        shared14 <- length(intersect(elementsList[[4]],elementsList[[1]]))
        shared124 <- length(intersect(intersect(elementsList[[1]],elementsList[[2]]),elementsList[[4]]))
        shared134 <- length(intersect(intersect(elementsList[[1]],elementsList[[3]]),elementsList[[4]]))
        shared234 <- length(intersect(intersect(elementsList[[2]],elementsList[[3]]),elementsList[[4]]))
        overlap1234 <- intersect(intersect(elementsList[[1]],elementsList[[2]]),intersect(elementsList[[3]],elementsList[[4]]))
        shared1234  <- length(overlap1234) 
    }
    if (ellipses ==5){
        area[5] <- length(elementsList[[5]])
        shared15 <- length(intersect(elementsList[[1]],elementsList[[5]]))
        shared25 <- length(intersect(elementsList[[2]],elementsList[[5]]))
        shared35 <- length(intersect(elementsList[[3]],elementsList[[5]]))
        shared45 <- length(intersect(elementsList[[4]],elementsList[[5]]))
        shared125 <- length(intersect(intersect(elementsList[[1]],elementsList[[2]]),elementsList[[5]]))
        shared135 <- length(intersect(intersect(elementsList[[1]],elementsList[[3]]),elementsList[[5]]))
        shared145 <- length(intersect(intersect(elementsList[[1]],elementsList[[4]]),elementsList[[5]]))
        shared235 <- length(intersect(intersect(elementsList[[2]],elementsList[[3]]),elementsList[[5]]))
        shared245 <- length(intersect(intersect(elementsList[[2]],elementsList[[4]]),elementsList[[5]]))
        shared345 <- length(intersect(intersect(elementsList[[3]],elementsList[[4]]),elementsList[[5]]))
        shared1235 <- length(intersect(intersect(elementsList[[1]],elementsList[[2]]),intersect(elementsList[[3]],elementsList[[5]]))) 
        shared1245 <- length(intersect(intersect(elementsList[[1]],elementsList[[2]]),intersect(elementsList[[4]],elementsList[[5]])))
        shared1345 <- length(intersect(intersect(elementsList[[1]],elementsList[[3]]),intersect(elementsList[[4]],elementsList[[5]])))
        shared2345 <- length(intersect(intersect(elementsList[[2]],elementsList[[3]]),intersect(elementsList[[4]],elementsList[[5]])))
        overlap12345 <- intersect(elementsList[[1]],intersect(intersect(elementsList[[2]],elementsList[[3]]),intersect(elementsList[[4]],elementsList[[5]])))
        shared12345  <-length(overlap12345)
      }
    }

  if (ellipses  == 2){
    venn.plot<-draw.pairwise.venn(length(elementsList[[1]]), length(elementsList[[2]]),shared12,
                                category = categories,scaled = scaleF,fill = colorValues, 
                                lty = c(rep("solid",2)), cex = intLabSize, cat.cex = 2.5,euler.d = euler)
    overlap <- overlap12

  }
  if (ellipses  == 3){
      venn.plot<-draw.triple.venn(area[1], area[2],area[3], 
                            shared12,shared23,shared13,shared123,
                            category = categories,scaled = scaleF,euler.d = euler,fill = colorValues, 
                            lty = c(rep("solid",3)), cex = intLabSize, cat.cex = 2.5)
      overlap <- overlap123
  }
  if (ellipses  == 4){
      venn.plot <- draw.quad.venn(area[1], area[2],area[3],area[4], 
                            shared12,shared13,shared14,shared23,shared24,
                            shared34,shared123,shared124,shared134,shared234,shared1234,
                            category = categories,scaled = scaleF,fill = colorValues, 
                            lty = c(rep("solid",4)), cex = intLabSize, cat.cex = 2.5,euler.d = euler)
      overlap <- overlap1234      
  }
  if (ellipses  == 5){
      venn.plot <- draw.quintuple.venn(area[1], area[2],area[3],area[4],area[5], 
                                 shared12,shared13,shared14,shared15,shared23,shared24,shared25,
                                 shared34,shared35,shared45,shared123,shared124,shared125,shared134,
                                 shared135,shared145,shared234,shared235,shared245,shared345,shared1234,
                                 shared1235,shared1245,shared1345,shared2345,shared12345,
                                 category = categories,scaled = scaleF,fill = colorValues, 
                                 lty = c(rep("solid",5)), cex = intLabSize, cat.cex = 2.5,euler.d = euler)
      overlap <- overlap12345
  }
  return(overlap)
}