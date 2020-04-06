# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("topGO")
# BiocManager::install("KEGGREST")
# BiocManager::install("org.Hs.eg.db")


library(topGO)
library(KEGGREST)
library(org.Hs.eg.db)
library(biomaRt)

# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html

getwd()

###################################################
# topGO 
###################################################

#######  Data preparation: ###########
# List of genes identifiers, gene scores, list of differentially expressed genes 
# or a criteriafor selecting genes based on their scores, 
# as well as gene-to-GO annotations are all collected and storedin a single R object.

## we found out that IC1 is related to sex
## IC9 is related to MLL cytogenetic subgroup
IC1 <- Jade_ICA1$S[,1]

# keep only genes with projection > 2.5
# genes_ic1 <- IC1[abs(IC1) > 2.5]
genes_ic1 <- IC1

# get gene symbols from TPM_resm file
ic1_symbols <- TPM_resm[TPM_resm$gene_id %in% names(genes_ic1), c('gene_id','Gene')]
ic1_symbols$proj <- genes_ic1[as.vector(ic1_symbols$gene_id)] 
  
# prepare geneList for topGO
geneList <- ic1_symbols$proj
names(geneList) <- ic1_symbols$Gene
# remove genes with no symbol name
geneList <- geneList[names(geneList) != ""]


# gene seletion function used in creating topGOData object
# keep only the genes with projection >2.5 in each IC
geneSelFunc <- function (score) {
  return(abs(score) > 2.5)
}

# Create topGOData object
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = geneSelFunc,
              nodeSize = 5,
              annot = annFUN.org, mapping = "org.Hs.eg.db",
              ID = "symbol")

# retrieve genes2GO list from the "expanded" annotation in GOdata
allGO = genesInTerm(GOdata)

####### Running the enrichment tests:##########
# Using the object created in the first step the user 
# can perform enrichmentanalysis using any feasible mixture 
# of statistical tests and methods that deal with the GO topology.

## Fisher’s exact test which is based on gene counts, 
## Kolmogorov-Smirnov liketest which computes enrichment based on gene scores
## each gene has ascore (representing how differentially expressed a gene is)
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "fisher") 
# algorithm can be :  "classic", "elim", "weight01"
# statistic can be : "fisher", "ks"


tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
head(tab, 10)
dim(tab)



# Figure shows the the subgraph induced by the 5 most significant GOterms 
# as identified by the elim algorithm.  Significant nodes are represented 
# as rectangles.  The plotted graphis the upper induced graph generated 
# by these significant nodes.
par(cex = 0.18) # adjust size of text in the graph
showSigOfNodes(GOdata, score(resultKS), 
               firstSigNodes = 10, useInfo ='all')


####### Gene Location on chromasome (useful for sex) ################
# check locations of the genes with high projections
genes_highproj <- genes_ic1[geneSelFunc(genes_ic1)] 
genes_highproj[names(sort(abs(genes_highproj),decreasing = TRUE))]
locs <- TPM_resm[TPM_resm$gene_id %in% names(genes_highproj), c("gene_id","Gene","Location")]
# add projection 
locs$proj <- genes_highproj[match(locs$gene_id, names(genes_highproj))]
# 10 genes with the highest projection
locs[order(abs(locs$proj), decreasing = TRUE),][1:10,]



# chromosome count for locs
chr_cnt <- table(substr(locs$Location,1,4))
chr_cnt

# check for HOXA gene in IC9 
# Several gene expression studies have found a characteristic transcriptional signature in MLL-F AMLs,
# which includes high expression of most HOXA genes
hoxa <- substr(locs$Gene,1,4) == "HOXA"
locs[hoxa,]

