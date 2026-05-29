##########################################################
# Summarization and Visualization of Gene-centric 
# Noncoding Analysis Results using STAARpipelineSummary
##########################################################

### Load required R packages
system("dx-restore-folder /STAARpipeline/R_packages/rstudio_workbench_ukbrap_trial.zilinli_iu.2025-01-26T02-24-06.tar.gz")

############## load source code
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(dplyr)
library(RColorBrewer)
palette("default")

###########################################################
#           User Input
###########################################################

## results path
input_path <- "/mnt/project/UKB_500K_WGS_staarpipeline/CardioAge/revised/Gene_Centric_NonCoding/"
output_path <- "/UKB_500K_WGS_staarpipeline/CardioAge/revised/STAARpipelineSummary/Gene_Centric_NonCoding/"

## alpha level
alpha <- 2.5E-06 #0.05/20000

## results summary
cMAC_cutoff <- 20

###########################################################
#           Main Function 
###########################################################
col_c <- rainbow(7)

## SPA status
# When there is an imbalance in case-control ratios for binary traits, Saddlepoint Approximations (SPA) may be used to obtain accurate p-values. 
use_SPA <- FALSE

results_noncoding_genome <- get(load(paste0(input_path,"results_noncoding_genome_batch.Rdata")))

#######################################################
#     summarize unconditional analysis results
#######################################################

# 7 masks
results_UTR_genome <- c()
results_upstream_genome <- c()
results_downstream_genome <- c()
results_promoter_CAGE_genome <- c()
results_promoter_DHS_genome <- c()
results_enhancer_CAGE_genome <- c()
results_enhancer_DHS_genome <- c()

for(kk in 1:length(results_noncoding_genome))
{
  print(kk)
  for(k in 1:7)
  {
    results <- results_noncoding_genome[[kk]][k][[1]]
    if(is.null(results)==FALSE)
    {
      ### UTR
      if(results[3]=="UTR")
      {
        results_UTR_genome <- rbind(results_UTR_genome,results)
      }
      
      ### upstream
      if(results[3]=="upstream")
      {
        results_upstream_genome <- rbind(results_upstream_genome,results)
      }
      
      ### downstream
      if(results[3]=="downstream")
      {
        results_downstream_genome <- rbind(results_downstream_genome,results)
      }
      
      ### promoter_CAGE
      if(results[3]=="promoter_CAGE")
      {
        results_promoter_CAGE_genome <- rbind(results_promoter_CAGE_genome,results)
      }
      
      ### promoter_DHS
      if(results[3]=="promoter_DHS")
      {
        results_promoter_DHS_genome <- rbind(results_promoter_DHS_genome,results)
      }
      
      ### enhancer_CAGE
      if(results[3]=="enhancer_CAGE")
      {
        results_enhancer_CAGE_genome <- rbind(results_enhancer_CAGE_genome,results)
      }
      
      ### enhancer_DHS
      if(results[3]=="enhancer_DHS")
      {
        results_enhancer_DHS_genome <- rbind(results_enhancer_DHS_genome,results)
      }
      
    }
  }
}


###### cMAC_cutoff
# UTR
results_UTR_genome <- results_UTR_genome[results_UTR_genome[,"cMAC"]>cMAC_cutoff,]
# upstream
results_upstream_genome <- results_upstream_genome[results_upstream_genome[,"cMAC"]>cMAC_cutoff,]
# downstream
results_downstream_genome <- results_downstream_genome[results_downstream_genome[,"cMAC"]>cMAC_cutoff,]
# promoter_CAGE
results_promoter_CAGE_genome <- results_promoter_CAGE_genome[results_promoter_CAGE_genome[,"cMAC"]>cMAC_cutoff,]
# promoter_DHS
results_promoter_DHS_genome <- results_promoter_DHS_genome[results_promoter_DHS_genome[,"cMAC"]>cMAC_cutoff,]
# enhancer_CAGE
results_enhancer_CAGE_genome <- results_enhancer_CAGE_genome[results_enhancer_CAGE_genome[,"cMAC"]>cMAC_cutoff,]
# enhancer_DHS
results_enhancer_DHS_genome <- results_enhancer_DHS_genome[results_enhancer_DHS_genome[,"cMAC"]>cMAC_cutoff,]


###### whole-genome results
# UTR
save(results_UTR_genome,file="UTR.Rdata")
system(paste0("dx upload UTR.Rdata --path ",output_path,"UTR.Rdata"))
# upstream
save(results_upstream_genome,file="upstream.Rdata")
system(paste0("dx upload upstream.Rdata --path ",output_path,"upstream.Rdata"))
# downstream
save(results_downstream_genome,file="downstream.Rdata")
system(paste0("dx upload downstream.Rdata --path ",output_path,"downstream.Rdata"))
# promoter CAGE
save(results_promoter_CAGE_genome,file="promoter_CAGE.Rdata")
system(paste0("dx upload promoter_CAGE.Rdata --path ",output_path,"promoter_CAGE.Rdata"))
# promoter DHS
save(results_promoter_DHS_genome,file="promoter_DHS.Rdata")
system(paste0("dx upload promoter_DHS.Rdata --path ",output_path,"promoter_DHS.Rdata"))
# enhancer CAGE
save(results_enhancer_CAGE_genome,file="enhancer_CAGE.Rdata")
system(paste0("dx upload enhancer_CAGE.Rdata --path ",output_path,"enhancer_CAGE.Rdata"))
# enahncer DHS
save(results_enhancer_DHS_genome,file="enhancer_DHS.Rdata")
system(paste0("dx upload enhancer_DHS.Rdata --path ",output_path,"enhancer_DHS.Rdata"))

if(use_SPA)
{
  ### UTR
  UTR_sig <- results_UTR_genome[results_UTR_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(UTR_sig,file="UTR_sig.csv")
  system(paste0("dx upload UTR_sig.csv --path ",output_path,"UTR_sig.csv"))
  
  noncoding_sig <- c()
  noncoding_sig <- rbind(noncoding_sig,UTR_sig)
  
  ### upstream
  upstream_sig <- results_upstream_genome[results_upstream_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(upstream_sig,file="upstream_sig.csv")
  system(paste0("dx upload upstream_sig.csv --path ",output_path,"upstream_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,upstream_sig)
  
  ### downstream
  downstream_sig <- results_downstream_genome[results_downstream_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(downstream_sig,file="downstream_sig.csv")
  system(paste0("dx upload downstream_sig.csv --path ",output_path,"downstream_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,downstream_sig)
  
  ### promoter_CAGE
  promoter_CAGE_sig <- results_promoter_CAGE_genome[results_promoter_CAGE_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(promoter_CAGE_sig,file="promoter_CAGE_sig.csv")
  system(paste0("dx upload promoter_CAGE_sig.csv --path ",output_path,"promoter_CAGE_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,promoter_CAGE_sig)
  
  ### promoter_DHS
  promoter_DHS_sig <- results_promoter_DHS_genome[results_promoter_DHS_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(promoter_DHS_sig,file="promoter_DHS_sig.csv")
  system(paste0("dx upload promoter_DHS_sig.csv --path ",output_path,"promoter_DHS_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,promoter_DHS_sig)
  
  ### enhancer_CAGE
  enhancer_CAGE_sig <- results_enhancer_CAGE_genome[results_enhancer_CAGE_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(enhancer_CAGE_sig,file="enhancer_CAGE_sig.csv")
  system(paste0("dx upload enhancer_CAGE_sig.csv --path ",output_path,"enhancer_CAGE_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,enhancer_CAGE_sig)
  
  ### enhancer_DHS
  enhancer_DHS_sig <- results_enhancer_DHS_genome[results_enhancer_DHS_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(enhancer_DHS_sig,file="enhancer_DHS_sig.csv")
  system(paste0("dx upload enhancer_DHS_sig.csv --path ",output_path,"enhancer_DHS_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,enhancer_DHS_sig)
  
  write.csv(noncoding_sig,file=paste0(output_path,"noncoding_sig.csv"))
  system(paste0("dx upload noncoding_sig.csv --path ",output_path,"noncoding_sig.csv"))
}else
{
  ### UTR
  UTR_sig <- results_UTR_genome[results_UTR_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(UTR_sig,file="UTR_sig.csv")
  system(paste0("dx upload UTR_sig.csv --path ",output_path,"UTR_sig.csv"))
  
  noncoding_sig <- c()
  noncoding_sig <- rbind(noncoding_sig,UTR_sig)
  
  ### upstream
  upstream_sig <- results_upstream_genome[results_upstream_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(upstream_sig,file="upstream_sig.csv")
  system(paste0("dx upload upstream_sig.csv --path ",output_path,"upstream_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,upstream_sig)
  
  ### downstream
  downstream_sig <- results_downstream_genome[results_downstream_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(downstream_sig,file="downstream_sig.csv")
  system(paste0("dx upload downstream_sig.csv --path ",output_path,"downstream_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,downstream_sig)
  
  ### promoter_CAGE
  promoter_CAGE_sig <- results_promoter_CAGE_genome[results_promoter_CAGE_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(promoter_CAGE_sig,file="promoter_CAGE_sig.csv")
  system(paste0("dx upload promoter_CAGE_sig.csv --path ",output_path,"promoter_CAGE_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,promoter_CAGE_sig)
  
  ### promoter_DHS
  promoter_DHS_sig <- results_promoter_DHS_genome[results_promoter_DHS_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(promoter_DHS_sig,file="promoter_DHS_sig.csv")
  system(paste0("dx upload promoter_DHS_sig.csv --path ",output_path,"promoter_DHS_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,promoter_DHS_sig)
  
  ### enhancer_CAGE
  enhancer_CAGE_sig <- results_enhancer_CAGE_genome[results_enhancer_CAGE_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(enhancer_CAGE_sig,file="enhancer_CAGE_sig.csv")
  system(paste0("dx upload enhancer_CAGE_sig.csv --path ",output_path,"enhancer_CAGE_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,enhancer_CAGE_sig)
  
  ### enhancer_DHS
  enhancer_DHS_sig <- results_enhancer_DHS_genome[results_enhancer_DHS_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(enhancer_DHS_sig,file="enhancer_DHS_sig.csv")
  system(paste0("dx upload enhancer_DHS_sig.csv --path ",output_path,"enhancer_DHS_sig.csv"))
  
  noncoding_sig <- rbind(noncoding_sig,enhancer_DHS_sig)
  
  write.csv(noncoding_sig,file=paste0(output_path,"noncoding_sig.csv"))
  system(paste0("dx upload noncoding_sig.csv --path ",output_path,"noncoding_sig.csv"))
}

# More relaxed threshold
alpha_1E5 <- 1E-5
if(use_SPA)
{
  ### UTR
  UTR_sig_1E5 <- results_UTR_genome[results_UTR_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(UTR_sig_1E5,file="UTR_sig_1E5.csv")
  system(paste0("dx upload UTR_sig_1E5.csv --path ",output_path,"UTR_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- c()
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,UTR_sig_1E5)
  
  ### upstream
  upstream_sig_1E5 <- results_upstream_genome[results_upstream_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(upstream_sig_1E5,file="upstream_sig_1E5.csv")
  system(paste0("dx upload upstream_sig_1E5.csv --path ",output_path,"upstream_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,upstream_sig_1E5)
  
  ### downstream
  downstream_sig_1E5 <- results_downstream_genome[results_downstream_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(downstream_sig_1E5,file="downstream_sig_1E5.csv")
  system(paste0("dx upload downstream_sig_1E5.csv --path ",output_path,"downstream_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,downstream_sig_1E5)
  
  ### promoter_CAGE
  promoter_CAGE_sig_1E5 <- results_promoter_CAGE_genome[results_promoter_CAGE_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(promoter_CAGE_sig_1E5,file="promoter_CAGE_sig_1E5.csv")
  system(paste0("dx upload promoter_CAGE_sig_1E5.csv --path ",output_path,"promoter_CAGE_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,promoter_CAGE_sig_1E5)
  
  ### promoter_DHS
  promoter_DHS_sig_1E5 <- results_promoter_DHS_genome[results_promoter_DHS_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(promoter_DHS_sig_1E5,file="promoter_DHS_sig_1E5.csv")
  system(paste0("dx upload promoter_DHS_sig_1E5.csv --path ",output_path,"promoter_DHS_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,promoter_DHS_sig_1E5)
  
  ### enhancer_CAGE
  enhancer_CAGE_sig_1E5 <- results_enhancer_CAGE_genome[results_enhancer_CAGE_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(enhancer_CAGE_sig_1E5,file="enhancer_CAGE_sig_1E5.csv")
  system(paste0("dx upload enhancer_CAGE_sig_1E5.csv --path ",output_path,"enhancer_CAGE_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,enhancer_CAGE_sig_1E5)
  
  ### enhancer_DHS
  enhancer_DHS_sig_1E5 <- results_enhancer_DHS_genome[results_enhancer_DHS_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(enhancer_DHS_sig_1E5,file="enhancer_DHS_sig_1E5.csv")
  system(paste0("dx upload enhancer_DHS_sig_1E5.csv --path ",output_path,"enhancer_DHS_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,enhancer_DHS_sig_1E5)
  
  write.csv(noncoding_sig_1E5,file=paste0(output_path,"noncoding_sig_1E5.csv"))
  system(paste0("dx upload noncoding_sig_1E5.csv --path ",output_path,"noncoding_sig_1E5.csv"))
}else
{
  ### UTR
  UTR_sig_1E5 <- results_UTR_genome[results_UTR_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(UTR_sig_1E5,file="UTR_sig_1E5.csv")
  system(paste0("dx upload UTR_sig_1E5.csv --path ",output_path,"UTR_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- c()
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,UTR_sig_1E5)
  
  ### upstream
  upstream_sig_1E5 <- results_upstream_genome[results_upstream_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(upstream_sig_1E5,file="upstream_sig_1E5.csv")
  system(paste0("dx upload upstream_sig_1E5.csv --path ",output_path,"upstream_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,upstream_sig_1E5)
  
  ### downstream
  downstream_sig_1E5 <- results_downstream_genome[results_downstream_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(downstream_sig_1E5,file="downstream_sig_1E5.csv")
  system(paste0("dx upload downstream_sig_1E5.csv --path ",output_path,"downstream_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,downstream_sig_1E5)
  
  ### promoter_CAGE
  promoter_CAGE_sig_1E5 <- results_promoter_CAGE_genome[results_promoter_CAGE_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(promoter_CAGE_sig_1E5,file="promoter_CAGE_sig_1E5.csv")
  system(paste0("dx upload promoter_CAGE_sig_1E5.csv --path ",output_path,"promoter_CAGE_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,promoter_CAGE_sig_1E5)
  
  ### promoter_DHS
  promoter_DHS_sig_1E5 <- results_promoter_DHS_genome[results_promoter_DHS_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(promoter_DHS_sig_1E5,file="promoter_DHS_sig_1E5.csv")
  system(paste0("dx upload promoter_DHS_sig_1E5.csv --path ",output_path,"promoter_DHS_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,promoter_DHS_sig_1E5)
  
  ### enhancer_CAGE
  enhancer_CAGE_sig_1E5 <- results_enhancer_CAGE_genome[results_enhancer_CAGE_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(enhancer_CAGE_sig_1E5,file="enhancer_CAGE_sig_1E5.csv")
  system(paste0("dx upload enhancer_CAGE_sig_1E5.csv --path ",output_path,"enhancer_CAGE_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,enhancer_CAGE_sig_1E5)
  
  ### enhancer_DHS
  enhancer_DHS_sig_1E5 <- results_enhancer_DHS_genome[results_enhancer_DHS_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(enhancer_DHS_sig_1E5,file="enhancer_DHS_sig_1E5.csv")
  system(paste0("dx upload enhancer_DHS_sig_1E5.csv --path ",output_path,"enhancer_DHS_sig_1E5.csv"))
  
  noncoding_sig_1E5 <- rbind(noncoding_sig_1E5,enhancer_DHS_sig_1E5)
  
  write.csv(noncoding_sig_1E5,file=paste0(output_path,"noncoding_sig_1E5.csv"))
  system(paste0("dx upload noncoding_sig_1E5.csv --path ",output_path,"noncoding_sig_1E5.csv"))
}

###########################################################################
#       Manhattan Plot
###########################################################################

############## noncoding
### UTR
results_STAAR <- results_UTR_genome[,c(1,2,dim(results_UTR_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

genes_info_manhattan <- dplyr::left_join(genes_info,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "UTR"

### upstream
results_STAAR <- results_upstream_genome[,c(1,2,dim(results_upstream_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "upstream"

### downstream
results_STAAR <- results_downstream_genome[,c(1,2,dim(results_downstream_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "downstream"

### promoter_CAGE
results_STAAR <- results_promoter_CAGE_genome[,c(1,2,dim(results_promoter_CAGE_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "promoter_CAGE"

### promoter_DHS
results_STAAR <- results_promoter_DHS_genome[,c(1,2,dim(results_promoter_DHS_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "promoter_DHS"

### enhancer_CAGE
results_STAAR <- results_enhancer_CAGE_genome[,c(1,2,dim(results_enhancer_CAGE_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "enhancer_CAGE"

### enhancer_DHS
results_STAAR <- results_enhancer_DHS_genome[,c(1,2,dim(results_enhancer_DHS_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "enhancer_DHS"

## ylim
noncoding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]])
min_y <- ceiling(-log10(noncoding_minp)) + 1

if(is.infinite(min_y)){
  genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]] <- apply(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]],2,function(x){
    x <- ifelse(x==0,1e-308,x)
    return(x)
  })
  min_y <- 308
}

pch <- c(0,1,2,3,4,5,6)

figure1 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$UTR,sig.level=alpha,pch=pch[1],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

figure2 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$upstream,sig.level=alpha,pch=pch[2],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

figure3 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$downstream,sig.level=alpha,pch=pch[3],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

figure4 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$promoter_CAGE,sig.level=alpha,pch=pch[4],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

figure5 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$promoter_DHS,sig.level=alpha,pch=pch[5],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

figure6 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$enhancer_CAGE,sig.level=alpha,pch=pch[6],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

figure7 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$enhancer_DHS,sig.level=alpha,pch=pch[7],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

print("Manhattan plot")

png("gene_centric_noncoding_manhattan.png", width = 9, height = 6, units = 'in', res = 600)

print(figure1)
print(figure2,newpage = FALSE)
print(figure3,newpage = FALSE)
print(figure4,newpage = FALSE)
print(figure5,newpage = FALSE)
print(figure6,newpage = FALSE)
print(figure7,newpage = FALSE)

dev.off()
system(paste0("dx upload gene_centric_noncoding_manhattan.png --path ",output_path,"gene_centric_noncoding_manhattan.png"))

##########################################################
#          Q-Q Plot
##########################################################

print("Q-Q plot")
cex_point <- 1

png("gene_centric_noncoding_qqplot.png", width = 8, height = 8, units = 'in', res = 600)

### UTR
## remove unconverged p-values
observed <- sort(genes_info_manhattan$UTR[genes_info_manhattan$UTR < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### upstream
## remove unconverged p-values
observed <- sort(genes_info_manhattan$upstream[genes_info_manhattan$upstream < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### downstream
## remove unconverged p-values
observed <- sort(genes_info_manhattan$downstream[genes_info_manhattan$downstream < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### promoter_CAGE
## remove unconverged p-values
observed <- sort(genes_info_manhattan$promoter_CAGE[genes_info_manhattan$promoter_CAGE < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### promoter_DHS
## remove unconverged p-values
observed <- sort(genes_info_manhattan$promoter_DHS[genes_info_manhattan$promoter_DHS < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### enhancer_CAGE
## remove unconverged p-values
observed <- sort(genes_info_manhattan$enhancer_CAGE[genes_info_manhattan$enhancer_CAGE < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=5, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### enhancer_DHS
## remove unconverged p-values
observed <- sort(genes_info_manhattan$enhancer_DHS[genes_info_manhattan$enhancer_DHS < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=6, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

legend("topleft",legend=c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),ncol=1,bty="o",box.lwd=1,pch=0:6,cex=1.5,text.font=2)

dev.off()
system(paste0("dx upload gene_centric_noncoding_qqplot.png --path ",output_path,"gene_centric_noncoding_qqplot.png"))

rm(list=setdiff(ls(), c("input_path", "output_path", "alpha", "cMAC_cutoff", "col_c", "use_SPA"))); gc()
print(ls())