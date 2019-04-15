library(devtools)
library(ggfortify)
library(randomForest)
library(rpart)
library(neuralnet)

# Dataset: pos/neg label, k-mer counts
d <- read.table("kmer_list.txt", header = TRUE, sep = '\t')
# Value matrix
d.vals <- d[,-1]
# Value matrix without all-0 columns
d.vals <- d.vals[, colSums(d.vals != 0) > 0]

# PCA with scaling
d.pca <- prcomp(d.vals, center = TRUE, scale. = TRUE)
autoplot(d.pca, data = d, colour = 'glycosite', 
         loadings = TRUE, loadings.label = TRUE)

# PCA WITHOUT scaling and centering
d.pca <- prcomp(d.vals)
autoplot(d.pca, data = d, colour = 'glycosite', 
         loadings = TRUE, loadings.label = TRUE)

# K-means
set.seed(1)
autoplot(kmeans(d.vals, 3), data = d)

# Distribution
d.pre_dist = as.matrix(colSums(d.vals > 0))
d.dist <- hist(d.pre_dist, plot = TRUE, breaks = 200, )
               #labels = colnames(d.vals))
plot(density(d.pre_dist))
write.csv(x, file="kmer_count.txt")

# Distribution POS
d.pos <- as.matrix(d[which(d$glycosite=='pos'), ])
# No 0 columns
d.pos <- d.pos[, colSums(d.pos != 0) > 0]
# Value matrix
d.posm <- d.pos[, -1]
d.posm <- as.matrix(colSums(d.pos > 0))
d.pdist <- hist(d.posm, plot = TRUE, breaks = 200)

# Distribution NEG
d.neg <- as.matrix(d[which(d$glycosite=='neg'), ])
# No 0 columns
d.neg <- d.neg[, colSums(d.neg != 0) > 0]
# Value matrix
d.negm <- d.neg[, -1]
d.negm <- as.matrix(colSums(d.neg > 0))
d.ndist <- hist(d.negm, plot = TRUE, breaks = 200)

set.seed(1000)
d.pos <- as.matrix(d[which(d$glycosite=='pos'), ])
d.neg <- as.matrix(d[which(d$glycosite=='neg'), ])
d.neg <- d.neg[1:2200, ]
subset <- rbind(d.neg, d.pos)
subset <- subset[sample(nrow(subset)), ]
# Select half for train
set.seed(1000)
train_idx <- sample(1:nrow(subset),2200,replace=FALSE)
train <- subset[train_idx,] # select all these rows
test <- subset[-train_idx,] # select all but these rows

# PCA with scaling
subset2 <- subset[, -1]
subset2 <- subset2[, colSums(subset2 != 0) > 0]
d.pca <- prcomp(as.numeric(subset2), center = TRUE, scale. = TRUE)
autoplot(d.pca, data = subset2, colour = 'glycosite')
         #loadings = TRUE, loadings.label = TRUE)

# Class tree
tree <- rpart(glycosite ~., data = data.frame(train), method = "class")
printcp(tree)
rsq.rpart(tree)
pred <- predict(tree, type = "class")
table(pred)

# Random Forest
d.rf <- randomForest(glycosite ~ ., data = subset)
varImpPlot(d.rf)
plot(d.rf)
d.rf$confusion
d.rf

# Random Forest train/test
d.rf2 <- randomForest(glycosite ~ ., data = train)
varImpPlot(d.rf2)
d.rf2$confusion
rf2t <- predict(d.rf2, test)
plot(d.rf2)
