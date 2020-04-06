# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

BiocManager::install("DOSE")


library('org.Hs.eg.db')
library(DOSE)

universe_symbol = names(geneList)
gene_symbol = as.character(locs$Gene)

# convert gene symbol to entrez id
universe <- mapIds(org.Hs.eg.db, universe_symbol, 'ENTREZID', 'SYMBOL')
gene <- mapIds(org.Hs.eg.db, gene_symbol, 'ENTREZID', 'SYMBOL')

# DO Enrichment Analysis
yy = enrichDO(gene, ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH",
         universe, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
         readable = FALSE)

summary(yy)

# enrichment analysis based on the Network of Cancer Genes database (http://ncg.kcl.ac.uk/)
xx <- enrichNCG(gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                universe,minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                readable = FALSE)
summary(xx)

