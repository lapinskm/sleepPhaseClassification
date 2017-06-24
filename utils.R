
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

library(simecol)#approxTime

reseampleData <- function(data, timestamps) {
  result <- approxTime(x = data, xout = timestamps, method = "constant", rule = 2)
}

readHypnogtamAndResample <- function(timestamps, fileName) {

  edfHeader   <- readEdfHeader (fileName)
  edfData     <- readEdfSignals(edfHeader)
  annotations <- edfData$annotations

  hypnogram_1 <- edfData$annotation[c("onset", "annotation")]

  # This file contains only samples at start of phase.
  # To resample them you have to have at least one sample at end of phase.
  # So we generate them.
  endsOfPhase <- hypnogram_1[-1, "onset"]-1

  hypnogram_2 <- hypnogram_1[1:length(endsOfPhase),]
  hypnogram_2$onset <- endsOfPhase

  result <- rbind(hypnogram_1,hypnogram_2)

  # Weh have to put it in order
  result <- result[with(result, order(onset, annotation)), ]

  names(result)    <- c("t","stage")
  rownames(result) <- NULL

  # Resampling needs numeric values
  stage_lvls   <- levels(as.factor(result$stage))
  result$stage <- as.numeric(as.factor(result$stage))

  result <- reseampleData(result, timestamps )
  # Go back to factorial values of stages
  result$stage <- as.factor(result$stage)
  result
}
