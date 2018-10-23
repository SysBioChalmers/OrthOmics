plotDistributions <- function(lcpm,x,lcpm2,x2){
  # With the following code we can make a plot showing the raw log2CPM reads, the dotted
  # line indicates zero logCPM (= 1 cpm)
  nsamples <- ncol(x)
  # Ignore warning, some samples will have identical color.
  col <- brewer.pal(nsamples, "Paired") 
  par(mfrow = c(1, 2))
  #Plot raw data
  plot(
    density(lcpm[, 1]),
    col = col[1],
    ylim = c(0, 0.3),
    las = 2,
    main = paste(organism, " Raw data ", length(lcpm[,1]), " RNA"),
    xlab = "Log2-cpm"
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
    ylim = c(0, 0.3),
    las = 2,
    main = paste(organism, " Filtered data", length(lcpm2[,1]), " RNA"),
    xlab = "Log2-cpm"
  )
  abline(v = 0, lty = 3)
  for (i in 2:nsamples) {
    den <- density(lcpm2[, i])
    lines(den$x, den$y, col = col[i])
  }
}