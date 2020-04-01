## use ICA to predict sex

# install.packages("e1071")
library(e1071)
library(caret)

# load in ICA data
# IC1 is highly correlated with sex
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

confusionMatrix(as.factor(pred), as.factor(y))
# Confusion Matrix and Statistics
# 
# Reference
# Prediction   0   1
# 0 233   0
# 1   1 329
# 
# Accuracy : 0.9982     
# 95% CI : (0.9901, 1)
# No Information Rate : 0.5844     
# P-Value [Acc > NIR] : <2e-16     
# 
# Kappa : 0.9963     
# 
# Mcnemar's Test P-Value : 1          
#                                      
#             Sensitivity : 0.9957     
#             Specificity : 1.0000     
#          Pos Pred Value : 1.0000     
#          Neg Pred Value : 0.9970     
#              Prevalence : 0.4156     
#          Detection Rate : 0.4139     
#    Detection Prevalence : 0.4139     
#       Balanced Accuracy : 0.9979     
#                                      
#        'Positive' Class : 0   

# find the error prediction
compare <- as.data.frame(cbind(pred, y))
compare[compare$pred != compare$y,]

# Il y avait 07H060 et 10H029, tous deux M, 
# mais leur karyotype indique qu’ils ont perdu 
# le chr Y dans leurs cellules cancéreuses (45,X,-Y)

# 08H062 a un karyotype compliqué : 
# 47,XX,+21[15]/46,XX,i(21)(q10)[2]/46,XX[3] 

#07H060,10H029,08H062 are the samples in question
# X02H003     X02H009     X02H017 are ctrl samples (all female)
x_3 <- x[c(rownames(x)[3:5],"X07H060","X10H029","X08H062"),]
indx <- match(colnames(x), locs$gene_id)
colnames(x_3) <- locs$Gene[indx]

x_3t <-as.data.frame( t(x_3))
x_3t$loc <- locs$Location[indx]
x_3t
x_3t[order(x_3t$loc, decreasing = TRUE),]

###################### randomForest ###################
# # install.packages("randomForest")
# library(randomForest)
# 
# rf_classifier = randomForest(y~.,data =x, 
#                              ntree=100, mtry=10,
#                              nodesize = 5,
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
# compareRF[compareRF$predictForest != compareRF$y,]



