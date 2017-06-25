library(graphics)

source("utils.R")

#TODO: remove hard-codes

setwd("~/Source/sleepPhaseClassification")
psgFileName <- "SC4001E0-PSG.edf"
hypFileName <- "SC4001EC-Hypnogram.edf"

modelData  <- readAndPreprocessModelData(psgFileName, hypFileName)

img   <- modelData
img$t   <- NULL
img$stage <- NULL
image(z = data.matrix(img), col=rainbow(20))

plot(x=modelData[["t"]],y=modelData[["stage"]],type="l")

hypnogram_matrix   = data.matrix(as.numeric(modelData[["stage"]]))
hypnogram_timeline = data.matrix(as.numeric(modelData[["t"]]))

image (hypnogram_matrix,   col = rainbow(10),)

# preprocessing is finished. Now we can try build some models.

dataSize=nrow(modelData)


x <- modelData
x$stage <-NULL


#normalize <- function(x) { return ((x - mean(x)) / sqrt(sum(x^2)) )}
#x.norm=apply(x,2, normalize)
#x = x.norm

y <- modelData$stage


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
1-sum(ylearn_pred == ylearn)/length(ylearn_pred)
1-sum(ytest_pred == ytest  )/length(ytest_pred)

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

classInds  <- class.ind(neuralModelData$stage)
classNames <- paste0("P_",colnames(classInds))
colnames(classInds) <- paste0("P_",colnames(classInds))
neuralModelData$stage=NULL

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
neuralModel <- neuralnet(f, neuralModelData, linear.output=FALSE, hidden = c(50.50), threshold = 0.05, lifesign="full",lifesign.step=100)

#compare results
ylearn_pred <- neuralnet::compute(neuralModel, xlearn )$net.result
ylearn_pred <- apply(ylearn_pred, FUN=which.max, MARGIN = 1)
ylearn      <- apply(ylearn, FUN=which.max, MARGIN = 1)
1-sum(rep(1,length(ylearn_pred)) [ylearn_pred == ylearn])/length(ylearn_pred)

ytest_pred <- neuralnet::compute(neuralModel, xtest )$net.result
ytest_pred <- apply(ytest_pred, FUN=which.max, MARGIN = 1)
ytest      <- apply(ytest, FUN=which.max, MARGIN = 1)
1-sum(rep(1,length(ytest_pred)) [ytest_pred == ytest])/length(ytest_pred)
