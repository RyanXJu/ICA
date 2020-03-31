# use mixing matrix to do prediction
library(randomForest)  
library(e1071)  
library(caret) 

A <- Jade_ICA1$A
dim(A)


tissue <- datTraits$tissue

remove <- is.na(tissue)

# prepare training data
x <- A[!remove,]
y <- datTraits$tissue[!remove]


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
# - best performance: 0.1170551  

## create model with best parameters
svm_model_after_tune <- svm(x,y, kernel="radial", cost=10, gamma=0.1)

summary(svm_model_after_tune)

pred <- round(predict(svm_model_after_tune,x))
system.time(predict(svm_model_after_tune,x))
table(pred,y)



## best performance of each trait
for (trait in colnames(datTraits)){
  # empty samples
  remove <- is.na(datTraits[,trait])
  
  # prepare training data
  x <- A[!remove,]
  y <- datTraits[,trait][!remove]
  
  ################# tune model ##########################
  svm_tune <- tune(svm, train.x=x, train.y=y, 
                   kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(0,.1,.2,.3,.4,.5,1,2)))
  print(paste(trait, svm_tune$best.performance, sep = " "))
}
#
# [1] "tissue 0.122687342483278"
# [1] "sex 0.11857832578375"
# [1] "Overall_Survival_Time_days 1060887.08284484"
# [1] "LIC.frequency.absolute 1.71302870318418e-10"
# [1] "LSC17 91952.6250715572"
# [1] "mRNAsi 0.00851598279456273"


#