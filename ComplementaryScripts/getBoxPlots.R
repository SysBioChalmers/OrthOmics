getBoxPlots <- function(x,x2,titleStr,resultsPath,org,omic){ 
lcpm1     <- cpm(x, log = TRUE)
lcpm2     <- cpm(x2, log = TRUE)
nsamples  <- ncol(x)
colors    <- brewer.pal(nsamples, "Paired")
# We then visualize this with two graphs, before and after normalization:
par(mfrow = c(1, 2))
titleStr  <- paste(organism, '_',length(x2[,1]), ' ',omic,': Unnormalised',sep='')
boxplot(lcpm1,  las = 2, names = colnames(lcpm1), varwidth =TRUE, alpha = 0.2, col = colors)
title(main = titleStr, ylab = "Log-cpm")

titleStr  <- paste(organism, '_',length(x2[,1]), ' ',omic,': Normalised',sep='')
boxplot(lcpm2,  las = 2, names = colnames(lcpm1), varwidth =TRUE, alpha = 0.2, col = colors)
title(main = titleStr, ylab = "Log-cpm")
}
