## Classification and Regression Training
library(caret)
library(e1071)


index <- createDataPartition(iris$Species, p=0.75, list=FALSE)

# Subset training set with index
iris.training <- iris[index,]

# Subset test set with index
iris.test <- iris[-index,]

# Overview of algos supported by caret
names(getModelInfo())


# Train a model
model_knn <- train(iris.training[, 1:4], iris.training[, 5], method='knn')
model_cart <- train(iris.training[, 1:4], iris.training[, 5], method='rpart2')


##### Prediction

# Predict the labels of the test set
predictions<-predict.train(object=model_knn,iris.test[,1:4], type="raw")

# Evaluate the predictions
table(predictions)

# Confusion matrix 
confusionMatrix(predictions,iris.test[,5])


#####################MORE

# Train the model with preprocessing
model_knn <- train(iris.training[, 1:4], iris.training[, 5], method='knn', preProcess=c("center", "scale"))

# Predict values
predictions<-predict.train(object=model_knn,iris.test[,1:4], type="raw")

# Confusion matrix
confusionMatrix(predictions,iris.test[,5])