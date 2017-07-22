library(edfReader) # readEdfHeader, readEdfSignals
library(signal)    # specgram
library(simecol)   # approxTime

downloadDataFiles <- function(baseUrl, filenames, downloadLocation) {
  for (filename in filenames) {
    download.file(paste0(baseUrl,          filename),
                  paste0(downloadLocation, filename) );
  }
}

normalize <- function(x) {
  return ((x - mean(x)) / sqrt(sum(x^2)) )
}

prepareSpectrogram <- function(signalValue, windowTime, overlapFactor, signalFreq)
{
  #FFT params
  windowSize <- round(windowTime * signalFreq)
  window  <-  hanning(windowSize)
  overlap <- windowSize * overlapFactor
  fftInterval <- (windowSize-overlap)/signalFreq

  #perform FFT analysis of signal
  spectrogram <- specgram(signalValue, Fs = signalFreq, overlap = overlap, window = window)
  remove(signalValue, overlap, windowSize, window)
  #plot(spectrogram)

  #Get raw data from spectrogram
  signalSpectral <- abs(spectrogram$S)

  colnames(signalSpectral) <- spectrogram$t
  rownames(signalSpectral) <- spectrogram$f

  signalSpectral
}

cutSpectogramIntoFreqRanges <- function(data, rangeBorders) {
  #cutof over last border
  data <- data[as.numeric(rownames(data)) < rangeBorders[length(rangeBorders)] ,]
  dim(data)

  freqs <- as.numeric(rownames(data)) # frequencies of spectrogram
  times <- as.numeric(colnames(data))

  result <- data.frame(t=as.numeric(colnames(data)))

  #split data to ranges
  beginOfRange <- 0
  for(endOfRange in rangeBorders) {
    rangeName=paste0("r",beginOfRange , "_" , endOfRange, "hz")

    featureCollum <- colSums(data[freqs >beginOfRange &
                                    freqs<=endOfRange, ])

    result[[rangeName]] <- normalize(featureCollum)
    beginOfRange <- endOfRange
  }
  result
}

timestamps <- function (n, f) {
  interval=1/f;
   ret <- (1:n) * f
}

readSignalData <- function(fileName) {

  #Load edf file
  edfHeader <- readEdfHeader (fileName)
  edfData   <- readEdfSignals(edfHeader)
  #this script is memory usage optimized
  # We are interested in EEG Fpz-Cz signal only.
  # It should be enough to detect sleep phase.

  signalValue  <-  edfData$`EEG Fpz-Cz`$signal
  signalFreq   <-  edfData$`EEG Fpz-Cz`$sRate

  timeInterval <- 1/signalFreq
  dataSize <- length(signalValue)

  timeStamps <- 0:(dataSize-1)*timeInterval
  result <- data.frame(signalValue, timeStamps)
}

readSpectralData <- function(fileName) {
  #Load edf file
  edfHeader <- readEdfHeader (fileName)
  edfData   <- readEdfSignals(edfHeader)
  #this script is memory usage optimized
  # We are interested in EEG Fpz-Cz signal only.
  # It should be enough to detect sleep phase.
  signalValue  <-  edfData$`EEG Fpz-Cz`$signal
  signalFreq   <-  edfData$`EEG Fpz-Cz`$sRate
  remove(edfData)
  signalSpectral <- prepareSpectrogram(signalValue = signalValue,
                                       windowTime = 30,
                                       overlapFactor=1/2,
                                       signalFreq = signalFreq)

  rangeBordersFrequencies <- c( 0.3, 4, 6, 8, 10.5, 12.5, 14, 17, 20, 25)

  result <- cutSpectogramIntoFreqRanges(data = signalSpectral,
                                        rangeBorders = rangeBordersFrequencies)
}

readHypnogtam <- function(fileName) {

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
  result
}

reseampleHypnogramData <- function(data, timestamps) {
  # Resampling needs numeric values
  stage_lvls   <- levels(as.factor(data$stage))
  data$stage <- as.numeric(as.factor(data$stage))
  #Resample
  result <- approxTime(x = data, xout = timestamps, method = "constant", rule = 2)
  for(i in 1:length(stage_lvls)) {
     result[result==i] = stage_lvls[i]
  }
  # Go back to factorial values of stages
  result$stage <- factor(result$stage,stage_lvls)
  levels(result$stage) <- stage_lvls
  #Make it string again
  result$stage <- as.character(result$stage)
  #result$stage <- gsub(result$stage, "Sleep stage ", "")
  result
}

readAndPreprocessModelData <- function(psgFileName, hypFileName) {
  modelData     <- readSpectralData(psgFileName)
  hypnogramData <- readHypnogtam(hypFileName)
  hypnogramData <- reseampleHypnogramData(hypnogramData, modelData$t)
  modelData[["stage"]] <- hypnogramData$stage
  modelData
}

transformStringStagesToFactor <- function(stages) {
  stages  <- as.factor(stages)
  levels(stages)=gsub(" ","_",levels(stages))
  stages
}

prepareDataFromFileList <- function(fileList, dataDirectory) {
  i = 1
  retVal <- readAndPreprocessModelData (paste0(dataDirectory,filelist_1$psg[i]),
                                        paste0(dataDirectory,filelist_1$hypnogram[i]))

  for(i in (2:nrow(filelist_1)) ) {
    modelData <- readAndPreprocessModelData (paste0(dataDirectory,filelist_1$psg[i]),
                                             paste0(dataDirectory,filelist_1$hypnogram[i]))
    retVal <- rbind(retVal, modelData)
  }
  retVal <- retVal[retVal$stage != "Movement time" & retVal$stage != "Sleep stage ?",]

  retVal$stage <-transformStringStagesToFactor(retVal$stage)
  retVal
}

prepareTimeDomainData<- function(psgFile, hypFile) {
  signal <- readSignalData(psgFile)
  hypnogram <- readHypnogtam (hypFile)
  hypnogram <- reseampleHypnogramData (hypnogram, signal$timeStamps)
  ret <- data.frame(hypnogram$t, signal$signalValue, hypnogram$stage)
  ret
}

prepareTimeDomainDataFromFileList <- function(fileList, dataDirectory) {
  psgFiles <- paste0(dataDirectory, filelist$psg)
  hypFiles <- paste0(dataDirectory, filelist$hypnogram)
  retVal<-prepareTimeDomainData(psgFiles[1],hypFiles[1])
  for(i in (2:nrow(filelist_1)) ) {
    modelData <- read.csv(paste0(dataDirectory,filelist_1$psg[i]))
    retVal <- rbind(retVal, modelData)
  }

  retVal <- retVal[retVal$stage != "Movement time" & retVal$stage != "Sleep stage ?",]
  retVal$stage <- transformStringStagesToFactor(retVal$stage)
  retVal;
}