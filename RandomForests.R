library(MASS)
library(tidymodels)
library(tidyr)
library(rpart)
library(randomForest)
library(foreach)
library(ggplot2)
library(geometry)
library(mlbench)

#function to train the model and get the error
get.error <- function(i){
  #train model
  bio.rf <- randomForest(class ~ ., type = 'classification', data = train_data, 
                         importance=TRUE, ntree = nt, mtry = features, replace = T)
  
  #extract the prediction for the train test
  y.hat <- predict(bio.rf, newdata = subset(test_data, select = -c(class)))
  
  #confusion matrix
  cm <- table(test_data$class, y.hat)
  #errors committed
  as.numeric((cm[2,1] + cm[1,2])/sum(cm))
}

#open dataset
data("biopsy")
df <- biopsy[-c(1)]
df <- na.omit(df) #remove missing values
df$class <- factor(ifelse(df$class == 'benign', "0", "1")) #transform in label


# Split the data into training and testing sets
data_split <- initial_split(df, prop = 0.75)
train_data <- training(data_split)
test_data <- testing(data_split)



####################################################################
############## RANDOM FOREST #######################################
# for each experiments, each 'combination run' is repeated 100 times,
# and the results are averaged

#random forest using random input (attributes) selection
nt <- 100 # number of trees in each round
nf <- c(1:(ncol(df)-1)) # number of features considered in each split
errors_feat <- c()

for(features in nf){
  result <- t(sapply(1:2, get.error))
  errors_feat[features] <- mean(result)
}

errors_feat
ggplot() + geom_point(aes(x = nf, y = errors_feat))



# one tree for each random forest
nt <- 1
#different number of features
nf <- c(1:(ncol(df)-1)) # number of features extracted in each tree
errors_feat <- c()

for(features in nf){
  result <- t(sapply(1:100, get.error))
  errors_feat[features] <- mean(result)
}

errors_feat
ggplot() + geom_point(aes(x = nf, y = errors_feat))


#linear combination of inputs
L <- 3 # number of variavles to cobine
r <- choose(ncol(df)-1,L)

#extract all the linear combinations from L features
lc <- data <- data.frame(matrix(NA,    # Create empty data frame
                                nrow = nrow(df),
                                ncol = 0))

combinations <- combn(colnames(subset(df, select = -c(class))), L, simplify = F)

for(i in combinations){
  coef <- c(runif(L, -1,1))
  
  v <- unlist(i)
  
  lc[, paste("LC",paste(v, collapse = ''), sep = '')] <- df[v] %>%
    apply(., 1, function(x) dot(x, coef, d = T))
  
}

lc$class <- df$class #labels for the new dataset

#random forest on the new dataset
# Split the data into training and testing sets
data_split <- initial_split(lc, prop = 0.75)
train_data <- training(data_split)
test_data <- testing(data_split)

nt <- 100
#different number of features
nf <- c(1:(ncol(df)-1)) # number of features extracted in each tree
errors_feat <- c()

for(features in nf){
  result <- t(sapply(1:100, get.error))
  errors_feat[features] <- mean(result)
}

errors_feat
ggplot() + geom_point(aes(x = nf, y = errors_feat))
