library(e1071)
data(iris)
attach(iris)
library(MASS)
library(caTools)

SVM_model <- svm(train_data, train_label, probability = TRUE)
pred_prob <- predict(SVM_model, test_data, decision.values = TRUE, probability = TRUE)
heatmap(attr(pred_prob, "probabilities"), Colv = NA, col=viridis(100))

library(party)
library(mlr)
ct <- "group1"
mat.test <- as.data.frame(t(mat))
mat.test$celltype <- group[rownames(mat.test)]==ct
## Define the task
## Specify the type of analysis (e.g. classification) and provide data and response variable
task <- makeClassifTask(data = mat.test, target = "celltype")

#R_MAX_NUM_DLLS=200
method <- 'cforest.importance'
fv <- generateFilterValuesData(task, method = method)
results <- fv$data
results <- results[order(results[,3], decreasing=TRUE),]
barplot(results$cforest.importance)