###########################################################
# Fit STAAR null model

###########################################################

### Load required R packages
system("dx-restore-folder /STAARpipeline/R_packages/rstudio_workbench_ukbrap_trial.zilinli_iu.2025-01-26T02-24-06.tar.gz")

rm(list=ls())
gc()

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Matrix)
library(STAARpipeline)

## read phenotype and covariates
data <- read.csv("/mnt/project/UKB_500K_WGS_staarpipeline/CardioAge/revised/UKB_catboost_All_Data_Prediction_resid_corrected_Final.csv")
colnames(data)[1] <- "userId"
colnames(data)[2] <- "age"
colnames(data)[3] <- "sex"

# define age2
data$age2 <- (data$age)^2


##Load PCs 
pcs = read.table("/mnt/project/nullmodels/Data/PC_GRM/sGRM_and_PCs/output.pca.score") 
colnames(pcs)[1:2] <- c("userId","userID2")
colnames(pcs)[3:22] <- paste0("PC",1:20)

## merge phenotype with pcs
fullDat = merge(data,pcs,by="userId")

# overlap with WES samples
######## Get UKB WGS 500K sample ID
gId <- get(load("/mnt/project/nullmodels/Data/Phenotype/sample_id/UKB_500K_WGS_sampleid.RData"))
length(gId)

## Remove missing 
fullDatRed = fullDat[(fullDat$userId%in%gId)&!is.na(fullDat$CardioAG_corr)&!is.na(fullDat$PC1)&!is.na(fullDat$age)&!is.na(fullDat$sex),] 

# 21734   114
dim(fullDatRed)


########## rank normal transformation
##### CardioAG

# run regression
fullDatRed$CardioAG.resid <- resid(lm(CardioAG_corr ~ sex+age+age2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = fullDatRed))

# rank normal
fullDatRed$CardioAG.norm <- sd(fullDatRed$CardioAG_corr)*scale(qnorm((rank(fullDatRed$CardioAG.resid,na.last="keep")-0.5)/length(fullDatRed$CardioAG.resid)))

# save data
save(fullDatRed,file="WGS_CardioAge.20250912.Rdata")
system("dx upload WGS_CardioAge.20250912.Rdata --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/WGS_CardioAge.20250912.Rdata")

#########################################################
#       Fit null model
#########################################################
## load phenotype
phenotype <- fullDatRed

## load GRM
sgrm <- get(load("/mnt/project/nullmodels/Data/PC_GRM/sGRM_and_PCs/output.sparseGRM.sGRM.RData"))
sample_id <- unlist(lapply(strsplit(colnames(sgrm),"_"),`[[`,1))

colnames(sgrm) <- sample_id
rownames(sgrm) <- sample_id

##### CardioAG

a <- Sys.time()
### fit null model
obj.STAAR.UKB.CardioAge <- fit_nullmodel(CardioAG.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = phenotype, kins = sgrm, kins_cutoff = 0.022, id = "userId", use_sparse = TRUE,family = gaussian(link = "identity"), verbose=T)
b <- Sys.time()
b - a

save(obj.STAAR.UKB.CardioAge,file="WGS_CardioAge_obj.STAAR.UKB.20260125.Rdata")
system("dx upload WGS_CardioAge_obj.STAAR.UKB.20260125.Rdata --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/WGS_CardioAge_obj.STAAR.UKB.20260125.Rdata")