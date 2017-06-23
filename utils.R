
normalize <- function(x) {
  return ((x - mean(x)) / sqrt(sum(x^2)) )
}

cutDataIntoFreqRanges <- function(data, ranges) {
  #cutof upper
  data <- data[as.numeric(rownames(data)) < ranges[length(ranges)] ,]
  dim(data)

  freqs <- as.numeric(rownames(data)) # frequencies of spectrogram
  times <- as.numeric(colnames(data))

  result <- data.frame(t=as.numeric(colnames(data)))

  #split data to ranges
  beginOfRange <- 0
  for(endOfRange in ranges) {
    rangeName=paste0("r",beginOfRange , "_" , endOfRange, "hz")

    featureCollum <- colSums(data[freqs >beginOfRange &
                                    freqs<=endOfRange, ])
    result[[rangeName]] <- normalize(featureCollum)
    beginOfRange <- endOfRange
  }
  result
}
