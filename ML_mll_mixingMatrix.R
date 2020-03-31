# use mixing matrix to do prediction
library(randomForest)  
library(e1071)  
library(caret) 

A <- Jade_ICA1$A
dim(A)

cyto_group <- as.data.frame(cyto_group)

mll <- cyto_group$MLL_tras

remove <- is.na(mll)

# prepare training data
x <- A[!remove,]
colnames(x) <- paste("IC", c(1:40), sep = "")
y <-cyto_group$MLL_tras[!remove]


################# tune model ##########################
svm_tune <- tune(svm, train.x=x, train.y=as.numeric(y), 
                 kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(0,.1,.2,.3,.4,.5,1,2)))

summary(svm_tune)
print(svm_tune)
# Parameter tuning of ‘svm’:
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   cost gamma
# 10   0.1
# 
# - best performance: 0.1170551  

## create model with best parameters
svm_model_after_tune <- svm(x,y, kernel="radial",type = 'C-classification',
                            cost=10, gamma=0.1)

summary(svm_model_after_tune)

pred <- predict(svm_model_after_tune,x)
system.time(predict(svm_model_after_tune,x))
table(pred,y)

## best performance of each trait
for (grp in colnames(cyto_group)){
  # empty samples
  remove <- is.na(cyto_group[,grp])
  
  # prepare training data
  x <- A[!remove,]
  y <- cyto_group[,grp][!remove]
  
  # skip subgroup is total sample is too small (<10)
  if (sum(y) < 10){
    print(grp)
    next
  }
  
  ################# tune model ##########################
  svm_tune <- tune(svm, train.x=x, train.y=y, type = 'C-classification',
                   kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(0,.1,.2,.3,.4,.5,1,2)))
  
  print(paste(grp, svm_tune$best.performance, sep = " "))
}




######### Naive Bayes ################################


# split data into training and testing, p=0.8
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

## Categorical data only:
nb <- naiveBayes(y_training ~ ., data = as.data.frame(training))

predict(nb, testing)
predict(nb, test_nb, type = "raw")

pred <- predict(nb, testing)
table(pred, y_testing)


########### CV Naive Bayes use caret #######################
train_control <- trainControl(method = "cv", number = 10)
nbGrid <- expand.grid(usekernel = c(TRUE, FALSE),
                      fL = 0:5,
                      adjust = seq(0, 5, by = 1))

nb_tune = train(training, as.factor(y_training), 'nb',
                trControl=train_control,
                tuneGrid = nbGrid)
# warnings are not a big deal; one of you new data points falls outside
# of the range of one of the distributions (probably the conditional).
# There's not much to do about that.
# https://github.com/topepo/caret/issues/793

nb_tune$results %>% top_n(5, wt = Accuracy) %>%arrange(desc(Accuracy))

# usekernel fL adjust  Accuracy     Kappa AccuracySD   KappaSD
# 1      TRUE  0      2 0.9478223 0.6613187 0.03085938 0.1781955
# 2      TRUE  1      2 0.9478223 0.6613187 0.03085938 0.1781955
# 3      TRUE  2      2 0.9478223 0.6613187 0.03085938 0.1781955
# 4      TRUE  3      2 0.9478223 0.6613187 0.03085938 0.1781955
# 5      TRUE  4      2 0.9478223 0.6613187 0.03085938 0.1781955
# 6      TRUE  5      2 0.9478223 0.6613187 0.03085938 0.1781955
plot(nb_tune)

# results for best model
confusionMatrix(nb_tune)

pred_nbcv <- predict(nb_tune, newdata = testing)
confusionMatrix(pred_nbcv, as.factor(y_testing))

# Confusion Matrix and Statistics
# 
# Reference
# Prediction FALSE TRUE
# FALSE   125    4
# TRUE      3    5
# 
# Accuracy : 0.9489          
# 95% CI : (0.8976, 0.9792)
# No Information Rate : 0.9343          
# P-Value [Acc > NIR] : 0.3159          
# 
# Kappa : 0.5611          
# 
# Mcnemar's Test P-Value : 1.0000          
#                                           
#             Sensitivity : 0.9766          
#             Specificity : 0.5556          
#          Pos Pred Value : 0.9690          
#          Neg Pred Value : 0.6250          
#              Prevalence : 0.9343          
#          Detection Rate : 0.9124          
#    Detection Prevalence : 0.9416          
#       Balanced Accuracy : 0.7661          
#                                           
#        'Positive' Class : FALSE    


###############  SVM CV Gridsearch  ###########
# https://dataaspirant.com/2017/01/19/support-vector-machine-classifier-implementation-r-caret-package/
numFolds <- trainControl(method = "cv", number = 10)
cGrid <- expand.grid( C=c(10^(-2:2)))
svmCV <- train(training, as.factor(y_training), 
               method = "svmLinear", trControl = numFolds, tuneGrid = cGrid)

svmCV
# Support Vector Machines with Linear Kernel 
# 
# 554 samples
# 40 predictor
# 2 classes: 'FALSE', 'TRUE' 
# 
# No pre-processing
# Resampling: Cross-Validated (10 fold) 
# Summary of sample sizes: 498, 499, 499, 499, 499, 498, ... 
# Resampling results across tuning parameters:
#   
#   C      Accuracy   Kappa    
# 1e-02  0.9837338  0.8558455
# 1e-01  0.9909740  0.9213302
# 1e+00  0.9782792  0.8396207
# 1e+01  0.9710065  0.7963545
# 1e+02  0.9710065  0.7963545
# 
# Accuracy was used to select the optimal model using the largest value.
# The final value used for the model was C = 0.1.
plot(svmCV)

prediction_svmCV <-predict(svmCV, newdata = testing)
table(y_testing, prediction_svmCV)

confusionMatrix(prediction_svmCV, as.factor(y_testing) )
# Confusion Matrix and Statistics
# 
# Reference
# Prediction FALSE TRUE
# FALSE   128    4
# TRUE      0    5
# 
# Accuracy : 0.9708         
# 95% CI : (0.9269, 0.992)
# No Information Rate : 0.9343         
# P-Value [Acc > NIR] : 0.04943        
# 
# Kappa : 0.7002         
# 
# Mcnemar's Test P-Value : 0.13361        
#                                          
#             Sensitivity : 1.0000         
#             Specificity : 0.5556         
#          Pos Pred Value : 0.9697         
#          Neg Pred Value : 1.0000         
#              Prevalence : 0.9343         
#          Detection Rate : 0.9343         
#    Detection Prevalence : 0.9635         
#       Balanced Accuracy : 0.7778         
#                                          
#        'Positive' Class : FALSE          
                                   
