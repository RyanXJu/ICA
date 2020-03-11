# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("MineICA")

library(MineICA)

set.seed(2004)
M <- matrix(rnorm(5000*6,sd=0.3),ncol=10)
M[1:100,1:3] <- M[1:100,1:3] + 2
M[1:200,1:3] <- M[1:200,4:6] +1

## Random initializations are used for each iteration of FastICA
## Estimates are clustered using hierarchical clustering with average linkage
res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
                          nbIt=3, funClus="hclust", method="average")

# Iq index measures the difference between the intra-clustersimilarity 
# and the extra-cluster similiarity.
res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
                          nbIt=3, funClus="hclust", method="ward")


Jade0 <- runICA(method = "JADE", X = M,nbComp = 2)
Jade0$A

Jade1 <- runICA(method = "JADE", X = M,nbComp = 2)
Jade1$A

fast0 <- runICA(method = "fastICA", X = M,nbComp = 2, alg.type = "deflation")
fast0$A

fast1 <- runICA(method = "fastICA", X = M,nbComp = 2, alg.type = "deflation")
fast1$A


###### buildIcaSetThis function builds an object of classIcaSet.
dat <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
rownames(dat) <- paste("g", 1:1000, sep="")
colnames(dat) <- paste("s", 1:10, sep="")

## build a data.frame containing sample annotations
annot <- data.frame(type=c(rep("a",5),rep("b",5)))
rownames(annot) <- colnames(dat)       

## run ICA
resJade <- runICA(X=dat, nbComp=3, method = "JADE")

## build params
params <- buildMineICAParams(resPath="toy/")

## build IcaSet object
icaSettoy <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
                         dat=dat, pData=annot, alreadyAnnot=TRUE)
params <- icaSettoy$params
icaSettoy <- icaSettoy$icaSet

# MainZ example
# BiocManager::install("breastCancerMAINZ")
library(breastCancerMAINZ)
data(mainz)
# how to use S4 object
# https://genomicsclass.github.io/book/pages/eset.html
mainz
exprs(mainz)[1:3,1:3]
pData(mainz)[1:3,4:7]

## run ICA
datExpr = exprs(mainz) # 22283*200
Sys.time()
resJade <- runICA(X=datExpr, nbComp=100, method = "JADE", maxit=10000)
Sys.time()

