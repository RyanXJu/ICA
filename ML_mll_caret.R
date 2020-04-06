## use SVM to predict MLL

# install.packages("e1071")
library(e1071)

# load in ICA data
# IC9 is related to MLL
IC1 <- Jade_ICA1$S[,9]
ensembl <- IC1[abs(IC1)>2.5]
ensembl

# get the expression of these genes in datExpr0
expr <- datExpr0[, names(ensembl)]
dim(expr)
expr[1:3,1:3]

mll <- cyto_group[,"MLL_tras"]

remove <- is.na(mll)

# prepare training data
x <- expr[!remove,]
y <- cyto_group[,"MLL_tras"][!remove]
y <- as.numeric(y)

# ############### create SVM model #######################
# svm_model <- svm(x,y)
# summary(svm_model)
# 
# pred <- round(predict(svm_model,x))
# # system.time(pred <- predict(svm_model,x))
# 
# table(pred,y)
# #           y
# # pred   0    1
# #   0   642   8
# #   1   1    40
# ################# tune model ##########################
# svm_tune <- tune(svm, train.x=x, train.y=y, kernel = "radial",
#                  ranges=list(cost=10^(-1:2), gamma=c(0,.1,.2,.3,.4,.5,1,2)))
# # kernel="radial" may over fit severly
# 
# summary(svm_tune)
# 
# # Parameter tuning of ‘svm’:
# #   
# #   - sampling method: 10-fold cross validation 
# # 
# # - best parameters:
# #   cost gamma
# # 10   0.1
# # 
# # - best performance: 0.06532576 
# 
# ## create model with best parameters
# svm_model_after_tune <- svm(x,y, cost=10, gamma=0.1)
# summary(svm_model_after_tune)
# 
# pred <- round(predict(svm_model_after_tune,x))
# system.time(predict(svm_model_after_tune,x))
# table(pred,y)
# #          y
# # pred   0    1
# #   0   643   0
# #   1   0     48
# 
# # find the error prediction
# compare <- as.data.frame(cbind(pred, y))
# compare[compare$pred != compare$y,]
# 
# 
# 
# ###################### randomForest ###################
# # install.packages("randomForest")
# library(randomForest)
# 
# rf_classifier = randomForest(y~.,data =x, 
#                              ntree=200, mtry=10,
#                              nodesize = 2,
#                              importance=TRUE)
# class(rf_classifier)
# str(rf_classifier)
# # importance of the features (genes)
# importance(rf_classifier)
# 
# predictForest <- round(predict(rf_classifier, newdata = x))
# table(y, predictForest)
# 
# compareRF <- as.data.frame(cbind(predictForest, y))
# compareRF[compareRF$predictForest != compare$y,]


#################### cross validation ####################
# install.packages("caret")
library(caret)


## data spliting ##
# https://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/

# split data into training and testing, p=0.75
set.seed(100)
inTrain <- createDataPartition(
  y = y,
  ## the outcome data are needed
  p = .8,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
str(inTrain)


training <- x[ inTrain,]
testing  <- x[-inTrain,]

y_training <- y[ inTrain]
y_testing <- y[-inTrain]

nrow(training)
nrow(testing)
length(y_training)
length(y_testing)

# # train rf model with training data
# rf_1 = randomForest(y_training~.,data =training, 
#                     method="class",
#                     ntree=200,
#                     nodesize = 5,
#                     importance=TRUE)
# class(rf_1)
# str(rf_1)
# # importance of the features (genes)
# importance(rf_1)
# 
# # performance
# print(rf_1)
# 
# pred_rf1<- round(predict(rf_1, newdata = testing))
# table(y_testing, pred_rf1)
# # pred_rf1
# # y_testing   0   1
# # 0 146   1
# # 1  16   9

## RF CV Gridsearch ###
numFolds <- trainControl(method = "cv", number = 10)
mtryGrid <- expand.grid(.mtry = c(10:sqrt(ncol(x)))) 

# mtryGrid <- expand.grid(.mtry = sqrt(ncol(x))) 

## weired problem, for train():
# the last column of the data used in train function need to be diagnosis.
cv_train <- cbind(training,y_training)

rfCV <- train(as.factor(y_training)~., data = cv_train, 
      method = "rf", trControl = numFolds, tuneGrid = mtryGrid)
rfCV
# Random Forest 
# 
# 553 samples
# 240 predictors
# 2 classes: '0', '1' 
# 
# No pre-processing
# Resampling: Cross-Validated (10 fold) 
# Summary of sample sizes: 498, 498, 497, 498, 498, 498, ... 
# Resampling results across tuning parameters:
#   
#   mtry  Accuracy   Kappa    
# 10    0.9800974  0.7795224
# 11    0.9819156  0.7870450
# 12    0.9801299  0.7979428
# 13    0.9819156  0.7870450
# 14    0.9819156  0.7870450
# 15    0.9819156  0.7870450
# 
# Accuracy was used to select the optimal model using the largest value.
# The final value used for the model was mtry = 11.
plot(rfCV)
prediction_rfCV <- predict(rfCV, newdata = testing)
table(y_testing, prediction_rfCV)
confusionMatrix(prediction_rfCV, as.factor(y_testing) )
# Confusion Matrix and Statistics
# 
# Reference
# Prediction   0   1
# 0 126   5
# 1   0   7
# 
# Accuracy : 0.9638          
# 95% CI : (0.9175, 0.9881)
# No Information Rate : 0.913           
# P-Value [Acc > NIR] : 0.01658         
# 
# Kappa : 0.7188          
# 
# Mcnemar's Test P-Value : 0.07364         
#                                           
#             Sensitivity : 1.0000          
#             Specificity : 0.5833          
#          Pos Pred Value : 0.9618          
#          Neg Pred Value : 1.0000          
#              Prevalence : 0.9130          
#          Detection Rate : 0.9130          
#    Detection Prevalence : 0.9493          
#       Balanced Accuracy : 0.7917          
#                                           
#        'Positive' Class : 0               
                                 

###############  SVM CV Gridsearch  ###########
# https://dataaspirant.com/2017/01/19/support-vector-machine-classifier-implementation-r-caret-package/
numFolds <- trainControl(method = "cv", number = 10)
cGrid <- expand.grid( C=c(10^(-2:2)))
svmCV <- train(as.factor(y_training)~., data = cv_train, 
      method = "svmLinear", trControl = numFolds, tuneGrid = cGrid)

svmCV
# 553 samples
# 240 predictors
# 2 classes: '0', '1' 
# 
# No pre-processing
# Resampling: Cross-Validated (10 fold) 
# Summary of sample sizes: 497, 497, 497, 497, 499, 499, ... 
# Resampling results across tuning parameters:
#   
#   C      Accuracy   Kappa    
# 1e-02  0.9908730  0.9003548
# 1e-01  0.9872354  0.8711368
# 1e+00  0.9872354  0.8711368
# 1e+01  0.9872354  0.8711368
# 1e+02  0.9872354  0.8711368
# 
# Accuracy was used to select the optimal model using the largest value.
# The final value used for the model was C = 0.01.

plot(svmCV)

prediction_svmCV <-predict(svmCV, newdata = testing)
table(y_testing, prediction_svmCV)

confusionMatrix(prediction_svmCV, as.factor(y_testing) )
# Confusion Matrix and Statistics
# 
# Reference
# Prediction   0   1
# 0 126   1
# 1   0  11
# 
# Accuracy : 0.9928          
# 95% CI : (0.9603, 0.9998)
# No Information Rate : 0.913           
# P-Value [Acc > NIR] : 4.993e-05       
# 
# Kappa : 0.9526          
# 
# Mcnemar's Test P-Value : 1               
#                                           
#             Sensitivity : 1.0000          
#             Specificity : 0.9167          
#          Pos Pred Value : 0.9921          
#          Neg Pred Value : 1.0000          
#              Prevalence : 0.9130          
#          Detection Rate : 0.9130          
#    Detection Prevalence : 0.9203          
#       Balanced Accuracy : 0.9583          
#                                           
#        'Positive' Class : 0      


### McNemar’s Test #########
# McNemar’s Test captures the errors made by both models. 
# Specifically, the No/Yes and Yes/No (A/B and B/A in your case) 
# cells in the confusion matrix. The test checks if there is a 
# significant difference between the counts in these two cells. That is all.
# Fail to Reject Null Hypothesis: Classifiers have a similar proportion of errors on the test set.
# Reject Null Hypothesis: Classifiers have a different proportion of errors on the test set.

# find the error prediction in testing data
names(prediction_svmCV) <- rownames(testing)
prediction_svmCV<-as.character(prediction_svmCV)

compare <- as.data.frame(cbind(prediction_svmCV, y_testing))
compare[compare$prediction_svmCV != compare$y_testing,]

#             prediction_svmCV y_testing
# X06H066                0         1
# allTraits[allTraits$sample_id == "X06H066",]

# find the error prediction in all data
prediction_svmCV_all <-predict(svmCV, newdata = x)
table(y, prediction_svmCV_all)

confusionMatrix(prediction_svmCV_all, as.factor(y) )
names(prediction_svmCV_all) <- rownames(x)
prediction_svmCV_all<-as.character(prediction_svmCV_all)

compare <- as.data.frame(cbind(prediction_svmCV_all, y))
compare[compare$prediction_svmCV_all != compare$y,]

# errors
# X02H033   ENAH
# X06H066   CASC5
# X11H095   
# X16H063

## all MLL samples
mll_samples<- allTraits[allTraits$cytogenetic.subgroup == "MLL translocations (+MLL FISH positive) (Irrespective of additional cytogenetic abnormalities)",
          c("sample_id", "status.at.sampling", "tissue",
            "karyotype", "MLL.PTD.mutation", "Transcriptome_Sequencer",
            "MLL.Partner", "RNASEQ_protocol")]

table(allTraits$MLL.Partner)

#        CASC5    ELL   ENAH    ENL   GAS7 MLLT10  MLLT3  MLLT4  MLLT6 SEPT_9 
# 654      1      3      1      4      1      5     11      8      1      2 

# X06H066 is the only one with CASC5
# 43~49,XX,del(3)(q12),+?6[2],del(8)(p11.2)[9],t(11;15)(q23;q14),-16[3],-20[3],+21[15],+22[2][cp20]