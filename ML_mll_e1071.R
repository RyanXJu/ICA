## use ICA to predict sex

# install.packages("e1071")
library(e1071)  # tune()
library(caret)  # train()

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

################# tune model ##########################
svm_tune <- tune(svm, train.x=x, train.y=y, 
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
# - best performance: 0.06529739 


## create model with best parameters
svm_model_after_tune <- svm(x,y, kernel="radial", cost=10, gamma=0.1)
# summary(svm_model_after_tune)
# Parameters:
#   SVM-Type:  eps-regression 
# SVM-Kernel:  radial 
# cost:  10 
# gamma:  0.1 
# epsilon:  0.1 

# Number of Support Vectors:  691
# *** high number of SV is not good, it is overfitting

pred <- round(predict(svm_model_after_tune,x))
system.time(predict(svm_model_after_tune,x))
table(pred,y)

confusionMatrix(as.factor(pred), as.factor(y))
# Confusion Matrix and Statistics
# 
# Reference
# Prediction   0   1
# 0 643   0
# 1   0  48
# 
# Accuracy : 1          
# 95% CI : (0.9947, 1)
# No Information Rate : 0.9305     
# P-Value [Acc > NIR] : < 2.2e-16  
# 
# Kappa : 1          
# 
# Mcnemar's Test P-Value : NA         
#                                      
#             Sensitivity : 1.0000     
#             Specificity : 1.0000     
#          Pos Pred Value : 1.0000     
#          Neg Pred Value : 1.0000     
#              Prevalence : 0.9305     
#          Detection Rate : 0.9305     
#    Detection Prevalence : 0.9305     
#       Balanced Accuracy : 1.0000     
#                                      
#        'Positive' Class : 0    

# find the error prediction
compare <- as.data.frame(cbind(pred, y))
compare[compare$pred != compare$y,]

