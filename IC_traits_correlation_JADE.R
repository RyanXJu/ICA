library(reshape2)
library(gplots)


setwd("/u/juxiao/AML_ICA")
getwd()


### load in fastICA results
jade <- load("jade_10000g40ICs_leucegene.RData")
jade

#########################
## clinical trait
#########################
traitData = read.csv("/u/juxiao/Datasets/leucegene_data/IRIC - Leucegene Database.csv",
                     stringsAsFactors=FALSE);
dim(traitData)
names(traitData)

traitData$sample_id <- paste("X",traitData$sample_id,sep = "")
allTraits <- traitData[traitData$sample_id %in% rownames(datExpr0),]

names(allTraits)
dim(allTraits)
###****** ADD LSC17 and mRNAsi data ******
fn_TPM_lsc17  = "/u/juxiao/LSC17/TPM_lsc17.tsv"
TPM_lsc17  = read.delim(fn_TPM_lsc17, as.is=TRUE, check.names=FALSE, header=FALSE) 
colnames(TPM_lsc17) <- c("sample_id", "LSC17")
TPM_lsc17$sample_id <- paste("X", TPM_lsc17$sample_id, sep="")

fn_TPM_mRNAsi  = "/u/juxiao/mRNAsi/TPM_StemScore.tsv"
TPM_mRNAsi = read.delim(fn_TPM_mRNAsi, as.is=TRUE, check.names=FALSE, header=FALSE)
colnames(TPM_mRNAsi) <- c("sample_id", "mRNAsi")
TPM_mRNAsi$sample_id <- paste("X", TPM_mRNAsi$sample_id, sep="")

allTraits <- merge(allTraits, TPM_lsc17, by="sample_id")
allTraits <- merge(allTraits, TPM_mRNAsi, by="sample_id")
dim(allTraits)

##################### prepare datTraits #########################
# select traits of intrest
# allTraits col: 18-FAB, 22-cytogenetic 30-OS_days 34-WHO
datTraits = allTraits[, c("tissue", "sex","Overall_Survival_Time_days",
                          "LIC.frequency.absolute","LSC17","mRNAsi") ]
rownames(datTraits) = allTraits[, 1]

datTraits$sex[datTraits$sex == "F"] <- 0
datTraits$sex[datTraits$sex == "M"] <- 1
datTraits$sex <- as.numeric(datTraits$sex)

datTraits$tissue[datTraits$tissue == "Bone marrow"] <- 0
datTraits$tissue[datTraits$tissue == "Blood" ] <- 1
datTraits$tissue <- as.numeric(datTraits$tissue)

dim(datTraits)
datTraits[1:3,1:6]

sample_sex <- !is.na(datTraits$sex)
sample_tissue <- !is.na(datTraits$tissue)

# # ICs vs sex
# ICA_sex_cor<-list()
# for (i in 1:50){
#   ICA_sex_cor<-c(ICA_sex_cor,
#                  cor(fast_ICA1$A[,i], datTraits$sex, use = "pairwise.complete.obs"))
# }
# ICA_sex_cor


# ICs vs traits
IC_traits_matrix <- c()
for (t in colnames(datTraits)) {
  IC_trait_cor<-c()
  for (i in 1:40){
    IC_trait_cor <-c(IC_trait_cor,
                     cor(Jade_ICA1$A[,i], datTraits[,t], use = "pairwise.complete.obs"))
  }
  IC_traits_matrix <- cbind(IC_traits_matrix, IC_trait_cor )
}

colnames(IC_traits_matrix) <- colnames(datTraits)
rownames(IC_traits_matrix) <- paste(rep("IC", 40) , c(1:40), sep="")
heatmap(IC_traits_matrix)
heatmap.2(IC_traits_matrix, col=redgreen(100), trace="none") 

###############################################################
##### cytogenetic
###############################################################
cytogenetic <- allTraits[,c("sample_id","cytogenetic.risk", "cytogenetic.subgroup")]
unique(cytogenetic.subgroup)
table(cytogenetic.subgroup)

cyto_group <- dcast(cytogenetic, sample_id ~ cytogenetic.subgroup)
rownames(cyto_group)<-cyto_group$sample_id
cyto_group<-cyto_group[2:19]
cyto_group <- !is.na(cyto_group)
# # rename cols
# [1] "MLL translocations (+MLL FISH positive) (Irrespective of additional cytogenetic abnormalities)"         
# [2] "t(15;17)(q24;q21)/PML-RARA (Irrespective of additional cytogenetic abnormalities)"                      
# [3] "Intermediate abnormal karyotype (except isolated trisomy/tetrasomy 8)"                                  
# [4] "Normal karyotype"                                                                                       
# [5] "Complex (3 and more chromosomal abnormalities)"                                                         
# [6] "Trisomy/tetrasomy 8 (isolated)"                                                                         
# [7] "Monosomy 5/ 5q-/Monosomy 7/ 7q- (less than 3 chromosomal abnormalities)"                                
# [8] "NUP98-NSD1(caryotype normal)"                                                                           
# [9] "t(8;21)(q22;q22)/RUNX1-RUNX1T1 (Irrespective of additional cytogenetic abnormalities)"                  
# [10] "inv(16)(p13.1q22)/t(16;16)(p13.1;q22)/CBFB-MYH11 (Irrespective of additional cytogenetic abnormalities)"
# [11] "EVI1 rearrangements (+EVI1 FISH positive) (Irrespective of additional cytogenetic abnormalities)"       
# [12] "t(6;9)(p23;q34) (Irrespective of additional cytogenetic abnormalities)"                                 
# [13] "germinal cell and hematological cancers with i(12)(p10)"                                                
# [14] "Monosomy17/del17p (less than 3 chromosomal abnormalities)"                                              
# [15] "Hyperdiploid numerical abnormalities only"                                                              
# [16] "Undetermined"                                                                                           
# [17] "NUP98 translocations (+NUP98 FISH positive) (Irrespective of additional cytogenetic abnormalities)"     
# [18] "germinal cell and hematological cancers with i(12)(p10)_not final" 

colnames(cyto_group) <- c("MLL_tras", "t15_17", "inter_abn",
                          "Normal", "Complex","Trisomy",
                          "Monosomy5", "NUP98_NSD1","t8_21",
                          "inv16", "EVI1", "t6_9",
                          "germinal", "Monosomy17", "Hyperdiploid",
                          "Undetermined", "NIP98","germinal_notfinal")

# ICs vs cytogroup
IC_cyto_matrix <- c()
for (c in colnames(cyto_group)) {
  IC_cyto_cor<-c()
  for (i in 1:40){
    IC_cyto_cor <-c(IC_cyto_cor,
                    cor(Jade_ICA1$A[,i], cyto_group[,c], use = "pairwise.complete.obs"))
  }
  IC_cyto_matrix <- cbind(IC_cyto_matrix, IC_cyto_cor )
}

colnames(IC_cyto_matrix) <- colnames(cyto_group)
rownames(IC_cyto_matrix) <- paste(rep("IC", 40) , c(1:40), sep="")
heatmap(IC_cyto_matrix)


heatmap.2(IC_cyto_matrix, col=redgreen(100), 
          density.info="none", trace="none") 




