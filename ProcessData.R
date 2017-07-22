library(graphics)
library(caret) #confusionMatrix

setwd("~/Source/sleepPhaseClassification")

source("utils.R")

baseUrl    <- readLines("link.txt")
filelist_1 <- read.csv("fileList_1.csv", stringsAsFactors = FALSE)
filelist_2 <- read.csv("fileList_2.csv", stringsAsFactors = FALSE)

dataDir    <- "data/"

#downloadDataFiles(baseUrl, filelist_1$psg,       dataDir)
#downloadDataFiles(baseUrl, filelist_1$hypnogram, dataDir)
#downloadDataFiles(baseUrl, filelist_2$psg,       dataDir)
#downloadDataFiles(baseUrl, filelist_2$hypnogram, dataDir)

modelData <- prepareDataFromFileList(filelist_1, dataDir)

img   <- modelData
img$t   <- NULL
img$stage <- NULL
image(z = data.matrix(img), col=rainbow(20))

plot(x=modelData[["t"]],y=modelData[["stage"]],type="l")

hypnogram_matrix   = data.matrix(as.numeric(modelData[["stage"]]))
hypnogram_timeline = data.matrix(as.numeric(modelData[["t"]]))

head(modelData[["stage"]])

image (hypnogram_matrix,   col = rainbow(6))
legend(grconvertX(0.5, "device"), grconvertY(1, "device"),
       levels(modelData$stage), fill = rainbow(6), xpd = NA)

# preprocessing is finished. Now we can try build some models.

dataSize <- nrow(modelData)


x <- modelData
x$stage <-NULL
unique(modelData$stage)

#normalize <- function(x) { return ((x - mean(x)) / sqrt(sum(x^2)) )}
#x.norm=apply(x,2, normalize)
#x = x.norm

y <- modelData$stage


#Split data to learning and test sets
testIdx  <- sample(dataSize, dataSize * 0.5)

xlearn <- x[-testIdx,]
xtest  <- x[ testIdx,]

ylearn <- y[-testIdx]
ytest  <- y[ testIdx]

#Some models needs to be feed by matrix
xlearn_matrix <- data.matrix(xlearn)
xtest_matrix  <- data.matrix(xtest )

ylearn_matrix <- data.matrix(ylearn)
ytest_matrix  <- data.matrix(ytest )

# try SVM model
library("e1071")
svm_model  <- svm(x=xlearn_matrix, y=ylearn)
svm_model
ytest_pred  <- predict(svm_model, xtest_matrix)
ylearn_pred <- predict(svm_model, xlearn_matrix)

confusionMatrix(ylearn_pred, ylearn)
confusionMatrix(ytest_pred, ytest)

#try random forest
library(randomForest)
forest <- randomForest(xlearn, ylearn, xtest, ytest, ntree = 300, keep.forest = TRUE)
forest

#try RBF model
library(nnet)      #class.inds

#some models needs boolean class indicators
classInds   <- class.ind(modelData$stage)
ylearn_inds <- classInds[-testIdx,]
ytest_inds  <- classInds[ testIdx,]

#try rbf neural network model
library(RSNNS)

ylearn_numeric_matrix <- as.numeric(as.factor(ylearn_matrix))
ytest_numeric_matrix <- as.numeric(as.factor(ytest_matrix))

rbfn.model <- RSNNS::rbf(x=xlearn_matrix,
                         y=ylearn_inds,
                         size  =160,   # number of centers, ie, number of neurons in hidden layer
                         maxit =1000, # max number of iterations to learn
                         linOut=TRUE) # linear activation function (otherwise logistic)

ylearn_pred <- predict(rbfn.model, xlearn_matrix)
ylearn_pred <- apply(ylearn_pred, FUN=which.max, MARGIN = 1)

ytest_pred <- predict(rbfn.model, xtest_matrix)
ytest_pred <- apply(ytest_pred, FUN=which.max, MARGIN = 1)

confusionMatrix(ylearn_pred, ylearn)
confusionMatrix(ytest_pred, ytest)


#try neuralnet model
library(neuralnet) #neuralnet

#neural network needs other format of data so we have to prepare and split once againg
neuralModelData   <- modelData
neuralModelData$t <- NULL #time is not used in this case

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

ytest_pred <- neuralnet::compute(neuralModel, xtest )$net.result
ytest_pred <- apply(ytest_pred, FUN=which.max, MARGIN = 1)
ytest      <- apply(ytest, FUN=which.max, MARGIN = 1)

confusionMatrix(ylearn_pred, ylearn)
confusionMatrix(ytest_pred,  ytest)
