plotDistributions <- function(x1,x2,measured,UB,colVals,reps){
  
  lcpm    <- cpm(x1,log=TRUE)
  lcpm2    <- cpm(x2,log=TRUE)

  # With the following code we can make a plot showing the raw log2CPM reads, the dotted
  # line indicates zero logCPM (= 1 cpm)
  nsamples <- ncol(lcpm)
  #Assign condition-specific color to samples
  col <- c()
  for (i in 1:length(reps)){
    col <- c(col,rep(colVals[i],reps[i]))
  }

  if (all(measured==' proteins')){
    xlabel <- 'counts'
  }else{xlabel <- 'Log2-cpm'}
  
  par(mfrow = c(1, 2))
  #Plot raw data
  plot(
    density(lcpm[, 1]),
    col = col[1],
    ylim = c(0, UB),
    las = 2,
    main = paste(organism, " Raw data ", length(lcpm[,1]), measured),
    xlab = xlabel
  )
  abline(v = 0, lty = 3)
  for (i in 2:nsamples) {
    den <- density(lcpm[, i])
    lines(den$x, den$y, col = col[i])
  }
  #Plot filtered data
  plot(
    density(lcpm2[, 1]),
    col = col[1],
    ylim = c(0, UB),
    las = 2,
    main = paste(organism, " Filtered data", length(lcpm2[,1]), measured),
    xlab = xlabel
  )
  abline(v = 0, lty = 3)
  for (i in 2:nsamples) {
    den <- density(lcpm2[, i])
    lines(den$x, den$y, col = col[i])
  }
}