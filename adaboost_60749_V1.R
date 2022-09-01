# Clear environment
########
rm(list=ls())
gc()
########

# Import libraries
########
library(dplyr)
library(ROCR)
library(ggplot2)
library(adabag)
########

# Import data
########
setwd("C:/Users/bhavi/Desktop/Research")

# Import data
class1 <- read.csv("GSE81861_CRC_tumor_epithelial_cells_FPKM.csv")
class2 <- read.csv("GSE81861_CRC_NM_epithelial_cells_FPKM.csv")
gene <- class1$X

# Transpose and change to data frame
class1 <- data.frame(t(data.matrix(class1[,-1])))
class2 <- data.frame(t(data.matrix(class2[,-1])))

colnames(class1) <- gene
colnames(class2) <- gene
########

# Data Pre-processing
########
# Shuffle rows
set.seed(100)
class1<- class1[sample(nrow(class1)),]
set.seed(100)
class2<- class2[sample(nrow(class2)),]

class1['y'] <- 1
class2['y'] <- 0
l1 <- nrow(class1)
l2 <- nrow(class2)

class_all <- rbind(class1,class2)

# Remove genes not expressed or equally expressed in all cells
removeZeroVar <- function(df){
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}
data <- removeZeroVar(class_all)

data['y'] <- class_all['y']

gc()
########

# Pre-allocate variables
########
df_genes <- data.frame(genes=colnames(data[,-ncol(data)]), beta = numeric(ncol(data)-1))
df_genes$beta <- 0
row.names(df_genes) <- df_genes$genes

k=1
j=l1+1
nfolds = 3
l1fold = floor(l1/nfolds)
l2fold = floor(l2/nfolds)

xgboost.val <- rep(0,nfolds)

# create parameter list
params <- list(
  eta = .1,
  max_depth = 5,
  min_child_weight = 2,
  subsample = .8,
  colsample_bytree = .9
)
#########

# K-Fold Cross validation starts here
########
for (i in 1:nfolds){
  
  
# Train Test Split
#######
  
  if(i==nfolds){
    test_rows <- c(k:l1,j:(l1+l2))
  } else {
    test_rows <- c(k:(k+l1fold-1),j:(j+l2fold-1))
  }
  
  
  test_data <- data[test_rows,]
  train_data <- data[-test_rows,]
  test_data$y <- as.factor(test_data$y)
  train_data$y <- as.factor(train_data$y)
  k=k+l1fold
  j=j+l2fold
#######


# AdaBoost
#######
start.time <- Sys.time()
set.seed(100)
adaboost.fit <- boosting( y~.
                         , data = train_data
                       , boos = TRUE
                       , mfinal = 10
)
end.time <- Sys.time()
comp.time <- end.time - start.time
#######

# AUC
#######
adaboost.probabilities <- predict(adaboost.fit, newdata = test_data)
adaboost.y_pred <- ifelse(adaboost.probabilities > 0.5, 1, 0)
adaboost.pred <- prediction(adaboost.y_pred, y_test)
adaboost.auc.perf <- performance(adaboost.pred, measure = "auc")
adaboost.val[i] <- adaboost.auc.perf@y.values
#######

}
#######