library(MineICA)
library(dplyr)
library(gplots)


setwd("/u/juxiao/AML_ICA")
getwd()

#########################
## gene expression data
#########################
# load rna-seq data (normalized)
fn_TPM_resm = "/u/juxiao/Datasets/leucegene_data/genes_TPM.unstranded.annotated.tsv"
TPM_resm = read.delim(fn_TPM_resm)
dim(TPM_resm)
# [1] 60609   697
# note: there are only 691 samples in the 

# Load the filterd gene info
genes<-load(file="/u/juxiao/AML_WGCNA/genes_kept.RData")
genes

# choose kept gene info to use
gene_keep <- keep_avelog 
length(gene_keep)

sample_id <- colnames(TPM_resm)[startsWith(colnames(TPM_resm), "X")]
datExpr0 <- as.data.frame(t(TPM_resm[,sample_id]))
names(datExpr0) <- TPM_resm$gene_id
datExpr0[1:3,1:3]
# ENSG00000000003.15 ENSG00000000005.6 ENSG00000000419.12
# X01H001               0.02              0.00              35.68
# X01H002               0.05              0.00              20.44
# X02H003               0.08              0.01              32.15

# keep only flitered genes, log2(x+1)
datExpr0 <- log2(datExpr0[,gene_keep]+1)
rownames(datExpr0) <- sample_id
dim(datExpr0)


################################
### keep 10k genes sorted by IQR
################################
IQR_genes <- apply(datExpr0, 2, IQR, na.rm = TRUE)

# sorting the genes by IQR, keep 10k genes with largest IQR
ICA_genes <- names(sort(IQR_genes, decreasing = TRUE)[1:10000])
ICA_genes
#** 9457 genes overlap if use filterByExpr **


dat_ICA <- t(datExpr0[,ICA_genes])
dim(dat_ICA)
dat_ICA[1:3,1:3]
# X01H001   X01H002  X02H003
# ENSG00000129824.16 0.21412481 8.5168425 0.000000
# ENSG00000229807.12 7.30405465 0.2986583 7.878603
# ENSG00000067048.17 0.05658353 5.7319978 0.000000

#########################
### JADE #### 

## https://www.rasch.org/rmt/rmt11b.htm
## "Lack of convergence is an indication that the data do not fit the model well, 
## because there are too many poorly fitting observations. 
## A data set showing lack of convergence can usually be rescued by setting aside 
## for separate study the person or item performances which contain these 
## unexpected responses.""
#########################
# X: rows:genes, cols:samples
# A: mixing matrix; rows: samples, cols: ICs
# S: source matrix; rows: gens, cols:ICs

#*********** JADE may not converge (try reduce number of genes, or nb of ICs)
Jade_ICA1 <- runICA(method = "JADE", X = dat_ICA ,nbComp = 40, tol=10^-6, maxit = 1000)
Jade_ICA2 <- runICA(method = "JADE", X = dat_ICA ,nbComp = 40, tol=10^-6, maxit = 1000)
fast_ICA10k40ics <- runICA(method = "fastICA", X = t(dat_ICA) ,nbComp = 40, tol=10^-6, maxit = 1000, 
                           alg.type =  "parallel",fun = "logcosh")

cor_ICA_Jade <-cor(Jade_ICA1$S, Jade_ICA2$S)

# JADE is parametric, ICs are identical in each replicates
heatmap(cor_ICA_Jade)
save(Jade_ICA1, Jade_ICA2, file = "jade_10000g40ICs_leucegene.RData")
jade <- load("jade_10000g40ICs_leucegene.RData")
jade

heatmap.2(cor(Jade_ICA1$S, fast_ICA10k40ics$S), col= redgreen(100))


###########################
## fastICA
##########################

# fastICA can handle more genes
# sorting the genes by IQR, keep 10k genes with largest IQR
ICA_genes <- names(sort(IQR_genes, decreasing = TRUE)[1:10000])
ICA_genes

dat_ICA <- t(datExpr0[,ICA_genes])
dim(dat_ICA)
dat_ICA[1:3,1:3]


## Random initializations are used for each iteration of FastICA
## Estimates are clustered using hierarchical clustering with average linkage
fast_ICA1 <- clusterFastICARuns(X=as.matrix(dat_ICA), nbComp=50, alg.type="deflation",
                                nbIt=3, funClus="hclust", method="average")

fast_ICA2 <- clusterFastICARuns(X=as.matrix(dat_ICA), nbComp=50, alg.type="deflation",
                                nbIt=3, funClus="hclust", method="average")

cor_ICA_fast <- cor(fast_ICA1$S, fast_ICA2$S)

heatmap(cor_ICA_fast)
heatmap.2(cor_ICA_fast, col=redgreen(100),trace="none")
# correlations matrix between ICs in two different fastICA runs.
# shows that, ICs found by fastICA are not identical, but exist
# a strong and unique association between them, 
# indicating the high reproducibility of the ICAresults 
# *** have to save the ICs to obtain same downstream analysis results


# 
# # Iq index measures the difference between the intra-clustersimilarity 
# # and the extra-cluster similiarity.
# res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
#                           nbIt=3, funClus="hclust", method="ward")
fast_ICA3 <- clusterFastICARuns(X=as.matrix(dat_ICA), nbComp=50, alg.type="deflation",
                          nbIt=3, funClus="hclust", method="ward")

fast_ICA4 <- clusterFastICARuns(X=as.matrix(dat_ICA), nbComp=50, alg.type="deflation",
                                nbIt=3, funClus="hclust", method="ward")

cor_ICA_fast_ward <- cor(fast_ICA3$S, fast_ICA4$S)
heatmap(cor_ICA_fast_ward)
heatmap.2(cor_ICA_fast_ward, col=redgreen(100))

# cor_ICA_fast_13 <- cor(fast_ICA1$S, fast_ICA3$S)
# heatmap(cor_ICA_fast_13)

save(fast_ICA1, fast_ICA2, fast_ICA3, fast_ICA4, file = "fastICA_10000g50ICs_leucegene.RData")
fast <- load("fastICA_10000g50ICs_leucegene.RData")
fast


### correlation between ICs found by fastICA and JADE
heatmap.2(cor(Jade_ICA1$S, fast_ICA1$S), col = redgreen(100))



## ------------- investigate IC1 of fastICA1 -------------
IC1 <- fast_ICA1$S[,1]
sort(abs(IC1),decreasing = TRUE)[1:20]
# most of the high projection genes in thic IC is related to sex

