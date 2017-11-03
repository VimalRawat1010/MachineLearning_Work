

ExpData <- read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/MA/eset.txt",
                      header= T, row.names=1)


filter <- as.logical(rowSums(ExpData) > 10)
ExpData.new <- ExpData[filter,]
ExpData.mat <- as.matrix(as.data.frame(lapply(ExpData.new, as.integer)))
ExpData.Labels <- read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/MA/MA_header.txt",
                             header= F, row.names=1)
ExpData.Labels <- factor(ExpData.Labels$V2)
levels(ExpData.Labels)
length(levels(ExpData.Labels))
group <- as.numeric(ExpData.Labels)
names(group) <- colnames(ExpData.new)

library(xgboost)  # the main algorithm
library(archdata) # for the sample dataset
library(caret)    # for the confusionmatrix() function (also needs e1071 package)
library(plyr)
library(dplyr)    # for some data preperation



mat.df <- as.data.frame(t(ExpData.mat))
mat.df$class <- group
colnames(mat.df) <- c(rownames(ExpData.new), "class")

set.seed(717)


#=================================================================
# Split Data
#=================================================================

# 70/30% split of the training data,
# Make split index
train_index <- sample(1:nrow(mat.df), nrow(mat.df)*0.70)
# Full data set
data_variables <- as.matrix(mat.df[,-1001])
data_label <- mat.df[,"class"]
data_matrix <- xgb.DMatrix(data = as.matrix(mat.df), label = data_label)
# split train data and make xgb.DMatrix
train_data   <- data_variables[train_index,]
train_label  <- data_label[train_index]
train_matrix <- xgb.DMatrix(data = train_data, label = train_label)
# split test data and make xgb.DMatrix
test_data  <- data_variables[-train_index,]
test_label <- data_label[-train_index]
test_matrix <- xgb.DMatrix(data = test_data, label = test_label)




numberOfClasses <- length(unique(mat.df$class)) +1
xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "booster" = "gbtree",
                   "num_class" = numberOfClasses)
nround    <- 20 # number of XGBoost rounds
cv.nfold  <- 10
# Fit cv.nfold * cv.nround XGB models and save OOF predictions

cv_model <- xgb.cv(params = xgb_params,
                   data = train_matrix, 
                   nrounds = nround,
                   nfold = cv.nfold,
                   verbose = FALSE,
                   prediction = TRUE)

#stopCluster(cl)


OOF_prediction <- data.frame(cv_model$pred) %>%
  mutate(max_prob = max.col(., ties.method = "last"),
         label = train_label)
head(OOF_prediction)

library(viridis)
library(gplots)
# confusion matrix
cf.mat <- confusionMatrix(factor(OOF_prediction$label), 
                factor(OOF_prediction$max_prob -1),
                mode = "everything")

heatmap.2(cf.mat$table, main=paste0("Accuracy: ", cf.mat$overall[1], "\n", "xgb.cv", " CV(10), round 10"), Colv=NA, Rowv=NA, scale ="row", col=viridis(100), cellnote=cf.mat$table, trace="none", notecol= "black")



bst_model <- xgb.train(params = xgb_params,
                       data = train_matrix,
                       nrounds = nround)

test_pred <- predict(bst_model, newdata = test_matrix)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses) %>% 
  t() %>% data.frame() %>% 
  mutate(label = test_label,
         max_prob = max.col(., "last"))
# confusion matrix of test set
confusionMatrix(factor(test_prediction$label),
                factor(test_prediction$max_prob -1),
                mode = "everything")

### Variable Importance
# get the feature real names
names <-  colnames(mat.df[,-21429])
# compute feature importance matrix
importance_matrix = xgb.importance(feature_names = names, model = bst_model)
head(importance_matrix)
gp = xgb.plot.importance(importance_matrix)
print(gp) 






























