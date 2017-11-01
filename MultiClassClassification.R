set.seed(0)
G <- 5
N <- 30
M <- 1000
initmean <- 5
initvar <- 20
mat <- matrix(rnorm(N*M*G, initmean, initvar), M, N*G)
rownames(mat) <- paste0('gene', 1:M)
colnames(mat) <- paste0('cell', 1:(N*G))
group <- factor(sapply(1:G, function(x) rep(paste0('group', x), N)))
names(group) <- colnames(mat)


library(viridis)

heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)



set.seed(0)
upreg <- 5
upregvar <- 10
ng <- 100

diff <- lapply(1:G, function(x) {
  diff <- rownames(mat)[(((x-1)*ng)+1):(((x-1)*ng)+ng)]
  mat[diff, group==paste0('group', x)] <<- mat[diff, group==paste0('group', x)] + rnorm(ng, upreg, upregvar)
  return(diff)
})
names(diff) <- paste0('group', 1:G)
heatmap(mat, Rowv=NA, Colv=NA, col=viridis(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)



diff2 <- lapply(2:(G-1), function(x) {
  y <- x+G
  diff <- rownames(mat)[(((y-1)*ng)+1):(((y-1)*ng)+ng)]
  mat[diff, group %in% paste0("group", 1:x)] <<- mat[diff, group %in% paste0("group", 1:x)] + rnorm(ng, upreg, upregvar)
  return(diff)
})

heatmap(mat, Rowv=NA, Colv=NA, col=viridis(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)



mat.df <- as.data.frame(t(mat))
mat.df$class <- rep(levels(group), each=30)




library("xgboost")  # the main algorithm
library("archdata") # for the sample dataset
library("caret")    # for the confusionmatrix() function (also needs e1071 package)
library("dplyr")    # for some data preperation
library(doSNOW)

mat.df$class <- as.factor(mat.df$class)
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



#=================================================================
# Train Model
#=================================================================

numberOfClasses <- length(unique(mat.df$class)) +1
xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "booster" = "gbtree",
                   "num_class" = numberOfClasses)
nround    <- 100 # number of XGBoost rounds
cv.nfold  <- 10
# Fit cv.nfold * cv.nround XGB models and save OOF predictions

#train_matrix <- subset(mat.train, select= -c(class))

cl <- makeCluster(4, type = "SOCK")

# Register cluster so that caret will know to train in parallel.
registerDoSNOW(cl)





cv_model <- xgb.cv(params = xgb_params,
                   data = train_matrix, 
                   nrounds = nround,
                   nfold = cv.nfold,
                   verbose = FALSE,
                   prediction = TRUE)

stopCluster(cl)


OOF_prediction <- data.frame(cv_model$pred) %>%
  mutate(max_prob = max.col(., ties.method = "last"),
         label = train_label)
head(OOF_prediction)


# confusion matrix
confusionMatrix(factor(OOF_prediction$label), 
                factor(paste0('group', OOF_prediction$max_prob -1)),
                mode = "everything")















