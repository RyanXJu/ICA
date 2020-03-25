## use SVM to predict sex

# install.packages("e1071")
library(e1071)

# load in ICA data
IC1 <- Jade_ICA1$S[,1]
ensembl <- IC1[abs(IC1)>2.5]
ensembl

# get the expression of these genes in datExpr0
expr <- datExpr0[, names(ensembl)]
dim(expr)
expr[1:3,1:3]


sex <- datTraits$sex

remove <- is.na(sex)

# prepare training data
x <- expr[!remove,]
y <- datTraits$sex[!remove]

############### create SVM model #######################
svm_model <- svm(x,y)
summary(svm_model)

pred <- predict(svm_model,x)
system.time(pred <- predict(svm_model,x))

table(pred,y)

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
# - best performance: 0.02083586 

## create model with best parameters
svm_model_after_tune <- svm(x,y, kernel="radial", cost=10, gamma=0.1)
summary(svm_model_after_tune)

pred <- round(predict(svm_model_after_tune,x))
system.time(predict(svm_model_after_tune,x))
table(pred,y)
#           y
# pred   0    1
# 0     233   0
# 1     1    329

# find the error prediction
compare <- as.data.frame(cbind(pred, y))
compare[compare$pred != compare$y,]


###################### randomForest ###################
# install.packages("randomForest")
library(randomForest)

rf_classifier = randomForest(y~.,data =x, ntree=100, mtry=2, importance=TRUE)
class(rf_classifier)
str(rf_classifier)

predictForest <- round(predict(rf_classifier, newdata = x))
table(y, predictForest)

compareRF <- as.data.frame(cbind(predictForest, y))
compareRF[compareRF$predictForest != compare$y,]
