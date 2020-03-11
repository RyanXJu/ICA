## https://www.biostars.org/p/102088/

library(biomaRt)
# select mart and data set
bm <- useMart("ensembl")
bm <- useDataset("mmusculus_gene_ensembl", mart=bm)

# Get ensembl gene ids and GO terms
EG2GO <- getBM(mart=bm, 
               attributes=c('ensembl_gene_id','hgnc_symbol','go_id'))

# examine result
head(EG2GO,15)

# Remove blank entries
EG2GO <- EG2GO[EG2GO$go_id != '',]

# convert from table format to list format
geneID2GO <- by(EG2GO$go_id,
                EG2GO$ensembl_gene_id,
                function(x) as.character(x))

# examine result
head(geneID2GO)

# terms can be accessed using gene ids in various ways
geneID2GO$ENSMUSG00000098488
# [1] "GO:0009395" "GO:0008152" "GO:0005829" "GO:0030659" "GO:0004620"
# [6] "GO:0004623" "GO:0046872" "GO:0005515"
geneID2GO[['ENSMUSG00000098488']]
# [1] "GO:0009395" "GO:0008152" "GO:0005829" "GO:0030659" "GO:0004620"
# [6] "GO:0004623" "GO:0046872" "GO:0005515"