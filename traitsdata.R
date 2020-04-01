library(reshape2)

setwd("/u/juxiao/AML_ICA/ICA")
getwd()

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

# sample_sex <- !is.na(datTraits$sex)
# sample_tissue <- !is.na(datTraits$tissue)


###############################################################
##### cytogenetic subgroup
###############################################################
cytogenetic <- allTraits[,c("sample_id","cytogenetic.risk", "cytogenetic.subgroup")]
unique(cytogenetic$cytogenetic.subgroup)

## rename subgroups
table(cytogenetic$cytogenetic.subgroup)
# Complex (3 and more chromosomal abnormalities) 
# 102 
cytogenetic[cytogenetic == "Complex (3 and more chromosomal abnormalities)"] <- "Complex"

# EVI1 rearrangements (+EVI1 FISH positive) (Irrespective of additional cytogenetic abnormalities) 
# 16 
cytogenetic[cytogenetic == "EVI1 rearrangements (+EVI1 FISH positive) (Irrespective of additional cytogenetic abnormalities)"] <- "EVI1"

# germinal cell and hematological cancers with i(12)(p10) 
# 1 
cytogenetic[cytogenetic == "germinal cell and hematological cancers with i(12)(p10)"] <- "germinal"

# germinal cell and hematological cancers with i(12)(p10)_not final 
# 1 
cytogenetic[cytogenetic == "germinal cell and hematological cancers with i(12)(p10)_not final"] <- "germinal_notfinal"

# Hyperdiploid numerical abnormalities only 
# 2 
cytogenetic[cytogenetic == "Hyperdiploid numerical abnormalities only"] <- "Hyperdiploid"

# Intermediate abnormal karyotype (except isolated trisomy/tetrasomy 8) 
# 97 
cytogenetic[cytogenetic == "Intermediate abnormal karyotype (except isolated trisomy/tetrasomy 8)"] <- "inter_abn"

# inv(16)(p13.1q22)/t(16;16)(p13.1;q22)/CBFB-MYH11 (Irrespective of additional cytogenetic abnormalities) 
# 32 
cytogenetic[cytogenetic == "inv(16)(p13.1q22)/t(16;16)(p13.1;q22)/CBFB-MYH11 (Irrespective of additional cytogenetic abnormalities)"] <- "inv16"

# MLL translocations (+MLL FISH positive) (Irrespective of additional cytogenetic abnormalities) 
# 48
cytogenetic[cytogenetic == "MLL translocations (+MLL FISH positive) (Irrespective of additional cytogenetic abnormalities)"] <- "MLL_tras"

# Monosomy 5/ 5q-/Monosomy 7/ 7q- (less than 3 chromosomal abnormalities) 
# 24 
cytogenetic[cytogenetic == "Monosomy 5/ 5q-/Monosomy 7/ 7q- (less than 3 chromosomal abnormalities)"] <- "Monosomy5"

# Monosomy17/del17p (less than 3 chromosomal abnormalities) 
# 2
cytogenetic[cytogenetic == "Monosomy17/del17p (less than 3 chromosomal abnormalities)"] <- "Monosomy17"

# Normal karyotype 
# 286 
cytogenetic[cytogenetic == "Normal karyotype"] <- "Normal"

# NUP98 translocations (+NUP98 FISH positive) (Irrespective of additional cytogenetic abnormalities) 
# 1 
cytogenetic[cytogenetic == "NUP98 translocations (+NUP98 FISH positive) (Irrespective of additional cytogenetic abnormalities)"] <- "NUP98_trans"

# NUP98-NSD1(caryotype normal) 
# 5 
cytogenetic[cytogenetic == "NUP98-NSD1(caryotype normal)"] <- "NUP98_NSD1"

# t(15;17)(q24;q21)/PML-RARA (Irrespective of additional cytogenetic abnormalities) 
# 30 
cytogenetic[cytogenetic == "t(15;17)(q24;q21)/PML-RARA (Irrespective of additional cytogenetic abnormalities)"] <- "t15_17"

# t(6;9)(p23;q34) (Irrespective of additional cytogenetic abnormalities) 
# 3 
cytogenetic[cytogenetic == "t(6;9)(p23;q34) (Irrespective of additional cytogenetic abnormalities)"] <- "t6_9"

# t(8;21)(q22;q22)/RUNX1-RUNX1T1 (Irrespective of additional cytogenetic abnormalities) 
# 20 
cytogenetic[cytogenetic == "t(8;21)(q22;q22)/RUNX1-RUNX1T1 (Irrespective of additional cytogenetic abnormalities)"] <- "t8_21"

# Trisomy/tetrasomy 8 (isolated) 
# 20 
cytogenetic[cytogenetic == "Trisomy/tetrasomy 8 (isolated)"] <- "Trisomy"

# Undetermined 
# 1 
cytogenetic[cytogenetic == "Undetermined"] <- "Undetermined"


cyto_group <- dcast(cytogenetic, sample_id ~ cytogenetic.subgroup)
rownames(cyto_group)<-cyto_group$sample_id
cyto_group<-cyto_group[2:19]
cyto_group <- !is.na(cyto_group)

colSums(cyto_group)


save(datTraits, cyto_group, file = "datTraits.RData")
