library(MineICA)
library(dplyr)

setwd("/u/juxiao/AML_ICA")
getwd()

#########################
## gene expression data
#########################
# load rna-seq data (normalized)
fn_TPM_resm = "/u/juxiao/Datasets/leucegene_data/genes_TPM.unstranded.annotated.tsv"
TPM_resm = read.delim(fn_TPM_resm)
# dim(TPM_resm)
# [1] 60609   697
# note: there are only 691 samples in the 

# Load the filterd gene info
genes<-load(file="/u/juxiao/AML_WGCNA/genes_kept.RData")
genes

# choose kept gene info to use
gene_keep_avelog <- keep_avelog 
length(gene_keep_avelog)

gene_keep_filterByExpr <- keep_filterByExpr
length(gene_keep_filterByExpr)

sample_id <- colnames(TPM_resm)[startsWith(colnames(TPM_resm), "X")]
datExpr0 <- as.data.frame(t(TPM_resm[,sample_id]))
names(datExpr0) <- TPM_resm$gene_id
dim(datExpr0)
datExpr0[1:3,1:3]
# ENSG00000000003.15 ENSG00000000005.6 ENSG00000000419.12
# X01H001               0.02              0.00              35.68
# X01H002               0.05              0.00              20.44
# X02H003               0.08              0.01              32.15

# keep only flitered genes, log2(x+1)
datExpr0_avelog <- log2(datExpr0[,gene_keep_avelog]+1)
rownames(datExpr0_avelog) <- sample_id
dim(datExpr0_avelog)

datExpr0_filterByExpr <- log2(datExpr0[,gene_keep_filterByExpr]+1)
rownames(datExpr0_filterByExpr) <- sample_id
dim(datExpr0_filterByExpr)

################################
### keep 10k genes sorted by IQR
################################
IQR_genes_avelog <- apply(datExpr0_avelog, 2, IQR, na.rm = TRUE)
# sorting the genes by IQR, keep 10k genes with largest IQR
ICA_genes_avelog <- names(sort(IQR_genes_avelog, decreasing = TRUE)[1:2000])
ICA_genes_avelog

dat_ICA_avelog <- t(datExpr0_avelog[,ICA_genes_avelog])
dim(dat_ICA_avelog)
dat_ICA_avelog[1:3,1:3]
# X01H001   X01H002  X02H003
# ENSG00000129824.16 0.21412481 8.5168425 0.000000
# ENSG00000229807.12 7.30405465 0.2986583 7.878603
# ENSG00000067048.17 0.05658353 5.7319978 0.000000

IQR_genes_filterByExpr <- apply(datExpr0_filterByExpr, 2, IQR, na.rm = TRUE)
# sorting the genes by IQR, keep 10k genes with largest IQR
ICA_genes_filterByExpr <- names(sort(IQR_genes_filterByExpr, decreasing = TRUE)[1:2000])
ICA_genes_filterByExpr

dat_ICA_filterByExpr <- t(datExpr0_filterByExpr[,ICA_genes_filterByExpr])
dim(dat_ICA_filterByExpr)
dat_ICA_filterByExpr[1:3,1:3]


#########################
### JADE # not converge with 2800 genes

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

Jade_ICA_avelog <- runICA(method = "JADE", X = dat_ICA_avelog ,nbComp = 50, tol=10^-6, maxit = 1000)
Jade_ICA_filterByExpr <- runICA(method = "JADE", X = dat_ICA_filterByExpr ,nbComp = 50, tol=10^-6, maxit = 1000)

cor_ICA_Jade <-cor(Jade_ICA_avelog$S, Jade_ICA_filterByExpr$S)

# JADE is parametric, ICs are identical in each replicates
heatmap(abs(cor_ICA_Jade))
