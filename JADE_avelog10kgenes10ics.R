library(MineICA)
library(dplyr)

library(gplots)
library(biomaRt)
library(xlsx)
options(java.parameters = "-Xmx8000m")

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
gene_keep <- keep_avelog 
length(gene_keep)

# ############  Convert Ensembl_gene_id to gene_symbol ################
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# ensembl <- useMart("ensembl")
# ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# listAttributes(ensembl) #check all the attributes
# gene_keep <- getBM(filters="ensembl_gene_id_version",
#                   attributes=c("ensembl_gene_id_version", "hgnc_symbol"),
#                   values=gene_keep,
#                   mart=mart)


sample_id <- colnames(TPM_resm)[startsWith(colnames(TPM_resm), "X")]
datExpr0 <- as.data.frame(t(TPM_resm[,sample_id]))
dim(datExpr0)
names(datExpr0) <- TPM_resm$gene_id
datExpr0[1:3,1:3]

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
dat_ICA <- t(datExpr0[,ICA_genes])
dim(dat_ICA)
dat_ICA[1:3,1:3]
# X01H001   X01H002  X02H003
# ENSG00000129824.16 0.21412481 8.5168425 0.000000
# ENSG00000229807.12 7.30405465 0.2986583 7.878603
# ENSG00000067048.17 0.05658353 5.7319978 0.000000

#*********** JADE may not converge (try reduce number of genes, or nb of ICs)
Jade_ICA1 <- runICA(method = "JADE", X = dat_ICA,
                    nbComp = 40, tol=10^-6, maxit = 1000)
Jade_ICA2 <- runICA(method = "JADE", X = dat_ICA,
                    nbComp = 40, tol=10^-6, maxit = 1000)
cor_ICA_Jade <-cor(Jade_ICA1$S, Jade_ICA2$S)

# JADE is parametric, ICs are identical in each replicates
heatmap(cor_ICA_Jade)
heatmap.2(cor_ICA_Jade, col= redgreen(100), trace = "none")

save(Jade_ICA1, Jade_ICA2, file = "jade_avelog10kgenes40ICs_leucegene.RData")
jade <- load("jade_avelog10kgenes40ICs_leucegene.RData")
jade


############################################################
### load clinical data
############################################################
traits <- load("datTraits.RData")
traits

###########################################################
### ICs vs traits
###########################################################

ICs_traits_cor <-function(ICA_A, traits){
  ### ICA_A is the mixing matrix (rows:samples, cols: ICs)
  ### traits is the clinical outcome matrix (rows: samples, cols: traits)
  n <- dim(ICA_A)[2]
  IC_traits_matrix <- c()
  for (t in colnames(traits)) {
    IC_trait_cor<-c()
    for (i in 1:n){
      IC_trait_cor <-c(IC_trait_cor,
                       cor(ICA_A[,i], traits[,t], use = "pairwise.complete.obs"))
    }
    IC_traits_matrix <- cbind(IC_traits_matrix, IC_trait_cor )
  }
  colnames(IC_traits_matrix) <- colnames(traits)
  rownames(IC_traits_matrix) <- paste(rep("IC", n) , c(1:n), sep="")
  return(IC_traits_matrix)
}

### ICs vs sex, tissue, OS, stemness
JADE_traits <- ICs_traits_cor(Jade_ICA1$A, datTraits)
heatmap(JADE_traits)
heatmap.2(JADE_traits, col=redgreen(100), trace="none") 
# # IC1 is highly correlated with sex
# sort(abs(Jade_ICA1$S[,1]), decreasing = TRUE)[1:30]


### ICs vs cytogenetic subgroup
JADE_cytogenetic <- ICs_traits_cor(Jade_ICA1$A, cyto_group)
heatmap(JADE_cytogenetic)
heatmap.2(JADE_cytogenetic, col=redgreen(100), trace="none") 


#########################################################################
### enriched genes in each IC 
#########################################################################

# choose a threshold to define which genes are enriched in each IC
# visualize threshold in histogram
hist(Jade_ICA1$S, breaks = 200)
abline(v=c(-2.5,2.5), col="red")

# check how many genes are considered enriched in each IC
colSums(abs(Jade_ICA1$S) > 2.5)

# # transform ensembl_id to gene symbol 
# # use BiomaRt
# genes_symbol <- getBM(filters="ensembl_gene_id_version",
#                     attributes=c("ensembl_gene_id_version","hgnc_symbol"),
#                     values=rownames(Jade_ICA1$S),
#                     mart=mart)


# # use Leucegene TPM_resm
genes_symbol <- TPM_resm[TPM_resm$gene_id %in% rownames(Jade_ICA1$S), c("gene_id", "Gene")]
colnames(genes_symbol) <- c("ensembl_gene_id_version", "hgnc_symbol")

for (i in 1:40){
  ic <- Jade_ICA1$S[,i]
  projection <- ic[abs(ic) >2.5]
  print(projection)
  
  ic_genes <- genes_symbol[genes_symbol$ensembl_gene_id_version %in% names(projection) ,]
  ic_genes$projection <- projection[c(as.vector(ic_genes$ensembl_gene_id_version))]
  
  write.xlsx(ic_genes,file='Jade_avelog10kgenes40ics_enrichedgenes1.xlsx',sheetName = paste("IC",i,sep = ""),append = T)
}


# xlsx::write.xlsx(file1,file='XXXXX.xlsx',sheetName = '1')
# xlsx::write.xlsx(file2,file='XXXXX.xlsx',sheetName = '2',append = T) 
