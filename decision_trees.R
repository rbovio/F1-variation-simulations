#Classification tree instructions - from https://www.guru99.com/r-decision-trees.html
onlypopcomparisons = read.csv("all_parameters.csv")
onlypopcomparisons2 = onlypopcomparisons

#######RE-EDIT below so that there is a new variable for simulations that didn't produce hybrids or only produced 1 hybrid

#remove simulations where no hybrids survived
onlypopcomparisons2 = onlypopcomparisons2[complete.cases(onlypopcomparisons2[ , 4]),]
#convert p2_cv, h_sd, and h_cv NA's to zero
for(i in 1:nrow(onlypopcomparisons2)){
  if(onlypopcomparisons2$p2_mean[i] == 0){
    onlypopcomparisons2$p2_cv[i] = 0
  }
  if(is.na(onlypopcomparisons2$h_sd[i]) == T & is.na(onlypopcomparisons2$h_cv[i]) == T){
    onlypopcomparisons2$h_sd[i] = 0
    onlypopcomparisons2$h_cv[i] = 0
  }
}

onlypopcomparisons2$cvcomparison = NA
for(i in 1:nrow(onlypopcomparisons2)){
  if(onlypopcomparisons2$h_cv[i] < onlypopcomparisons2$p1_cv[i] & onlypopcomparisons2$h_cv[i] < onlypopcomparisons2$p2_cv[i]){
    onlypopcomparisons2$cvcomparison[i] = 'hybridcvsmaller'
  }
  else{
    onlypopcomparisons2$cvcomparison[i] = 'hybridcvlarger'
  }
}


# #From newoutput, need to extract only the rows where stat=="cvcomparison"
# onlypopcomparisons <- newoutput[newoutput[, "stat"] == "cvcomparison", ]
#I should shuffle the data around, because it is sorted
shuffle_index <- sample(1:nrow(onlypopcomparisons2))
onlypopcomparisons2 <- onlypopcomparisons2[shuffle_index, ]
#Do descriptive stats - what percent of observations is the CV smaller?
nrow(onlypopcomparisons2[onlypopcomparisons2$cvcomparison=="hybridcvsmaller",])/nrow(onlypopcomparisons2)

#Install and attach packages
#install.packages("rpart")
library(rpart)
#install.packages("rpart.plot")
library(rpart.plot)
#Create training/testing datasets
create_train_test <- function(data, size = 0.8, train = TRUE) {
  n_row = nrow(data)
  total_row = size * n_row
  train_sample <- 1:total_row
  if (train == TRUE) {
    return (data[train_sample, ])
  } else {
    return (data[-train_sample, ])
  }
}
data_train <- create_train_test(data = onlypopcomparisons2, size = 0.8, train = TRUE)
data_test <- create_train_test(data = onlypopcomparisons2, size = 0.8, train = FALSE)
#Build the model
fit <- rpart(cvcomparison ~ linkage + gene_action + p1_af + p2_af + esize + ss + ns, data = data_train, method="class")
summary(fit)
#Plot the model
rpart.plot(fit, extra = 106, fallen.leaves = T)
#Do sanity checks. Does the proportion of CV smaller you calculated earlier match
#the proportion in the top node? If they match - that is good. If they don't, 
#something is messed up.

#Make Predictions
predict_unseen <- predict(fit, data_test, type = "class")
#Test hybrids that had smaller CV and larger CV than both parents
table_mat <- table(data_test$cvcomparison, predict_unseen)
table_mat
#Sanity check - the values in this table should add to the number of observations
#in data_test

#Measure performance with a confusion matrix
accuracy_Test <- sum(diag(table_mat))/sum(table_mat)
print(paste('Accuracy for test', accuracy_Test))

