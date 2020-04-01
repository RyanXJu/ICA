library(MineICA)
library(dplyr)
library(gplots)
library(biomaRt)
library(xlsx)
options(java.parameters = "-Xmx8000m")

setwd("/u/juxiao/AML_ICA")
getwd()


#########################
## load ICA result from leucegene_ICA.R
#########################
jade <- load("jade_avelog10kgenes40ICs_leucegene.RData")
jade


############################################################
### load clinical data
############################################################
traits <- load("datTraits.RData")
traits

# colSums(cyto_group)
# Complex              EVI1          germinal germinal_notfinal 
# 102                16                 1                 1 
# Hyperdiploid         inter_abn             inv16          MLL_tras 
# 2                97                32                48 
# Monosomy17         Monosomy5            Normal        NUP98_NSD1 
# 2                24               286                 5 
# NUP98_trans            t15_17              t6_9             t8_21 
# 1                30                 3                20 
# Trisomy      Undetermined 
# 20                 1 
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

# IC9 seems highly correlated with MLL
JADE_cytogenetic[, "MLL_tras"]

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
  
  write.xlsx(ic_genes,file='Jade_avelog10kgenes40ics_enrichedgenes.xlsx',sheetName = paste("IC",i,sep = ""),append = T)
}


# xlsx::write.xlsx(file1,file='XXXXX.xlsx',sheetName = '1')
# xlsx::write.xlsx(file2,file='XXXXX.xlsx',sheetName = '2',append = T) 
