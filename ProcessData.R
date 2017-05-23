library(edfReader)
library(graphics)
library(signal)

#TODO: remove hard-codes
setwd("Source/SleepPhaseClassification")
fileName <- "SC4001E0-PSG.edf"

#Load edf files
edfHeader <- readEdfHeader (fileName)
edfData   <- readEdfSignals(edfHeader)
#this script is memory usage optimized
remove(edfHeader,fileName)

# We are interested in EEG Fpz-Cz signal only.
# It should be enough to detect sleep phase.
signalValue  <-  edfData$`EEG Fpz-Cz`$signal
signalFreq   <-  edfData$`EEG Fpz-Cz`$sRate
remove(edfData)

#FFT params
windowSize <- 30 * signalFreq # [30 s]
window  <-  hanning(windowSize)
#overlap <- window_size-1 #imposible due huge memory usage (it would generate fft value for each 1/100s)
overlap <- windowSize/2 # is default one - means there is value per 15 sec
fftInterval <- (windowSize-overlap)/signalFreq

#perform FFT analysis of signal
spectrogram <- specgram(signalValue, Fs = signalFreq, overlap = overlap, window = window)
remove(signalValue, overlap, windowSize, window)
#plot(spectrogram)

#Get raw data from spectrogram
signalSpectral <- abs(spectrogram$S)
ncol(signalSpectral)
nrow(signalSpectral)

colnames(signalSpectral) <- spectrogram$t
rownames(signalSpectral) <- spectrogram$f
length(spectrogram$t)
length(spectrogram$f)

remove(spectrogram)

#remove data values over 25 Hz - These values are meaningless for our needs - deleting them saves more memory
upperCutoffFrequency <- 25
signalSpectral <- signalSpectral[as.numeric(rownames(signalSpectral)) < upperCutoffFrequency ,]


nrow(signalSpectral)
ncol(signalSpectral)

#empiricaly selected aggregation range borders. Sums of these ranges are our features.

rangeBordersFrequencies <- c( 0.3, 4, 6, 8, 10.5, 12.5, 14, 17, 20, upperCutoffFrequency)

#rangeBordersFrequencies <- c( 0.3, 2, 3, 4, 5, 6, 7, 8, 9, 10.5,10, 11.5,12, 12.5, 13, 14,15.5, 17, 20, upperCutoffFrequency)
remove(upperCutoffFrequency)

modelData <- data.frame(t=as.numeric(colnames(signalSpectral)))

#split data to ranges
beginOfRange <- 0
for(endOfRange in rangeBordersFrequencies) {
  rangeName=paste0("r",beginOfRange , "_" , endOfRange, "hz")

  featureCollum <- colSums(signalSpectral[as.numeric(rownames(signalSpectral))>beginOfRange &
                                          as.numeric(rownames(signalSpectral))<=endOfRange, ])
  modelData[[rangeName]] <- featureCollum
  beginOfRange <- endOfRange
}
remove(signalSpectral, beginOfRange, endOfRange, featureCollum, rangeName, rangeBordersFrequencies)

#load hypnogram - modelled class
#TODO: remove hard-codes
hypnogram <- read.csv("SC4001E0-PSG.edf.csv",header = TRUE, stringsAsFactors = FALSE)
hypnogram[["beginTime"]]=0
hypnogram[["endTime"]]=0

#Calculate time of phase's begin and end
for (i in 2:nrow(hypnogram)) {
  hypnogram[["endTime"]][i-1]=hypnogram[["beginTime"]][i-1]+hypnogram[["duration"]][i-1]
  hypnogram[["beginTime"]][i]=hypnogram[["endTime"]][i-1]
}
hypnogram[["endTime"]][i]=hypnogram[["beginTime"]][i]+hypnogram[["duration"]][i]

modelData[["phase"]]="e"

# assign phase to each data point
firstRecordOfNextPhase=1
for (i in 1:nrow(hypnogram)) {
  phase = hypnogram[i,]
  for(j in firstRecordOfNextPhase : nrow(modelData)) {
    if  (modelData[["t"]][j] >= phase[["beginTime"]]  &&
         modelData[["t"]][j] <= phase[["endTime"]]) {
        modelData[["phase"]][j] = phase[["Sleep_stage"]]
    }
    else {
      firstRecordOfNextPhase = j
      break
    }
  }
}
remove(hypnogram, phase, i, j, firstRecordOfNextPhase)

modelData[["phase"]] = as.factor(modelData[["phase"]])

head(modelData)
tail(modelData)

# preprocessing is finished. Now we can try build some models.

dataSize=nrow(modelData)


x <- modelData
x$phase <-NULL


#normalize <- function(x) { return ((x - mean(x)) / sqrt(sum(x^2)) )}
#x.norm=apply(x,2, normalize)
#x = x.norm

y <- modelData$phase


#Split data to learning and test sets
testIdx  <- sample(dataSize, dataSize * 0.1)

xlearn <- x[-testIdx,]
xtest  <- x[ testIdx,]

ylearn <- y[-testIdx]
ytest  <- y[ testIdx]

#Some models needs to be feed by matrix
xlearn_matrix <- data.matrix(xlearn)
xtest_matrix  <- data.matrix(xtest )

# try SVM model
library("e1071")
svm_model  <- svm(x=xlearn_matrix, y=ylearn)
svm_model
ytest_pred  <- predict(svm_model, xtest_matrix)
ylearn_pred <- predict(svm_model, xlearn_matrix)
# compare result of prediction for learn and test data
1-sum(rep(1,length(ylearn_pred)) [ylearn_pred == ylearn])/length(ylearn_pred)
1-sum(rep(1,length(ytest_pred)) [ytest_pred == ytest] )/length(ytest_pred)

#try random forest
library(randomForest)
forest <- randomForest(xlearn, ylearn, xtest, ytest, ntree = 600, keep.forest = TRUE)
forest


#try neuralnet model
library(neuralnet)
library(nnet)#for class.ind only

#neural network needs other format of data so we have to prepare and split once againg
neuralModelData   <- modelData
neuralModelData$t <- NULL #time is not used in this case

classInds  <- class.ind(neuralModelData$phase)
classNames <- paste0("P_",colnames(classInds))
colnames(classInds) <- paste0("P_",colnames(classInds))
neuralModelData$phase=NULL

xlearn <- neuralModelData[-testIdx,]
ylearn <- classInds[-testIdx,]
xtest <- neuralModelData[ testIdx,]
ytest <- classInds[ testIdx,]

neuralModelData <- cbind(neuralModelData, classInds)

#create formula
lf <- paste(colnames(ylearn), collapse='+')
rf <- paste(colnames(xlearn), collapse='+')
f <- as.formula(paste(lf,'~',rf))

#teach model
neuralModel <- neuralnet(f, neuralModelData, linear.output=FALSE, hidden = c(40), threshold = 0.05, lifesign="full",lifesign.step=100)

#compare results
ylearn_pred <- neuralnet::compute(neuralModel, xlearn )$net.result
ylearn_pred <- apply(ylearn_pred, FUN=which.max, MARGIN = 1)
ylearn      <- apply(ylearn, FUN=which.max, MARGIN = 1)
1-sum(rep(1,length(ylearn_pred)) [ylearn_pred == ylearn])/length(ylearn_pred)

ytest_pred <- neuralnet::compute(neuralModel, xtest )$net.result
ytest_pred <- apply(ytest_pred, FUN=which.max, MARGIN = 1)
ytest      <- apply(ytest, FUN=which.max, MARGIN = 1)
1-sum(rep(1,length(ytest_pred)) [ytest_pred == ytest])/length(ytest_pred)
