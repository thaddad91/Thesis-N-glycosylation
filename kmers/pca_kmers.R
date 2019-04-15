#library(ggfortify)
#library(tidyverse)
#library(rpart)
#library(randomForest)
#library(e1071)
library(caret)
#library(kernlab)
#library(verification)
# Read 3-mer list
filename <- "kmer_list.txt"
d <- read.table(filename, sep="\t", row.names = NULL, header=TRUE)
# Separate pos from neg
d.pos <- d[which(d$glycosite=='pos'), ]
d.neg <- d[which(d$glycosite=='neg'), ]
# Subset neg for balance with pos
d.neg_sub <- sample_n(d.neg, 2200)
# Combine and remove neg columns
d.sub <- rbind(d.pos, d.neg_sub)
d.sub <-  d.sub[, colSums(d.sub != 0) > 0]
# Train and test set
#set.seed(10)
#tr.pos <- sample_n(d.pos, 1600)
#tr.neg <- sample_n(d.neg_sub, 1600)
#train <- rbind(tr.pos, tr.neg)
#train <-  train[, colSums(train != 0) > 0]
#te.pos <- sample_n(d.pos, 600)
#te.neg <- sample_n(d.neg_sub, 600)
#test <- rbind(te.pos, te.neg)
#test <-  test[, colSums(test != 0) > 0]

set.seed(100)
# Cross-validation
ctrl <- trainControl(method = "cv", savePred=T, classProb=T, number=10, repeats=3)
mod <- train(glycosite ~ ., data = d.sub, method = "neuralnet", trControl = ctrl)
head(mod$pred)
mod

mod2 <- train(glycosite ~ ., data = d.sub, method = "knn", trControl = ctrl)
head(mod2$pred)
mod2

# PCA
#d.pca <- prcomp(d.sub[, -1], center = TRUE, scale. = TRUE)
#autoplot(d.pca, data = d.sub, colour = 'glycosite')
         #loadings = TRUE, loadings.label = TRUE)

# Classification tree
#d.tree <- rpart(glycosite ~ ., data = d.sub, method = "class")
#printcp(d.tree)
#plotcp(d.tree)
# Plot tree
#plot(d.tree)
#text(d.tree)

# RandomForest
#d.rf <- randomForest(glycosite ~ ., data = d.sub)

# SVM
#d.svm = svm(as.factor(glycosite) ~ ., data = d.sub, scale = TRUE, kernel = "polynomial", cost = 5)
#table(predict(d.svm), d.sub$glycosite, dnn=c("Prediction", "Actual"))

#d.svm2 = svm(as.factor(glycosite) ~ ., data = train, scale = TRUE, kernel = "polynomial", cost = 5)
#table(predict(d.svm2), train$glycosite, dnn=c("Prediction", "Actual"))
#table(predict(d.svm2, test), test$glycosite, dnn=c("Prediction", "Actual"))
#plot(d.svm, d.sub, as.factor(glycosite) ~ ., fill = TRUE, grid = 50, slice = list(),
#     symbolPalette = palette(), svSymbol = "x", dataSymbol = "o")
