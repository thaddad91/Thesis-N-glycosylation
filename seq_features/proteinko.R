library(ggfortify)
library(randomForest)
data <- read.csv("table_pI", header = FALSE, sep = "\t", dec = ".")
# Balanced pos/neg set 
pos <- data[which(data$V1=='pos'), ]
neg <- data[which(data$V1=='neg'), ]
neg <- data[sample(nrow(data), 2200),]
d <- rbind(pos, neg)

# PCA
pca <- prcomp(d[, -1], center = TRUE, scale. = FALSE)
autoplot(pca, data = d, colour = 'V1', x = 5, y = 6)  # Other PCA's
#loadings = TRUE, loadings.label = TRUE)

# Randomforest
rf <- randomForest(V1 ~ ., data = d)
