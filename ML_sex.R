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

# Il y avait 07H060 et 10H029, tous deux M, 
# mais leur karyotype indique qu’ils ont perdu 
# le chr Y dans leurs cellules cancéreuses (45,X,-Y)



# 08H062 a un karyotype compliqué : 
# 47,XX,+21[15]/46,XX,i(21)(q10)[2]/46,XX[3] 

#07H060,10H029,08H062
x_3 <- x[c(rownames(x)[3:5],"X07H060","X10H029","X08H062"),]
indx <- match(colnames(x), locs$gene_id)
colnames(x_3) <- locs$Gene[indx]

x_3t <-as.data.frame( t(x_3))
x_3t$loc <- locs$Location[indx]
x_3t

###################### randomForest ###################
# install.packages("randomForest")
library(randomForest)

rf_classifier = randomForest(y~.,data =x, 
                             ntree=100, mtry=10,
                             nodesize = 5,
                             importance=TRUE)
class(rf_classifier)
str(rf_classifier)
# importance of the features (genes)
importance(rf_classifier)

predictForest <- round(predict(rf_classifier, newdata = x))
table(y, predictForest)

compareRF <- as.data.frame(cbind(predictForest, y))
compareRF[compareRF$predictForest != compare$y,]


#################### cross validation ####################
# install.packages("caret")
library(caret)
# https://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/

# split data into training and testing, p=0.75
set.seed(100)
inTrain <- createDataPartition(
  y = y,
  ## the outcome data are needed
  p = .75,
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

# train rf model with training data
rf_1 = randomForest(y_training~.,data =training, 
                             ntree=100, mtry=10,
                             nodesize = 5,
                             importance=TRUE)
class(rf_1)
str(rf_1)
# importance of the features (genes)
importance(rf_1)

# performance
print(rf_1)

pred_rf1<- round(predict(rf_1, newdata = testing))
table(y_testing, pred_rf1)

