##########################################################
# Summarization and Visualization of Gene-centric 
# Coding Analysis Results using STAARpipelineSummary
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
input_path <- "/mnt/project/UKB_500K_WGS_staarpipeline/CardioAge/revised/Gene_Centric_Coding_Incl_PTV/"
output_path <- "/UKB_500K_WGS_staarpipeline/CardioAge/revised/STAARpipelineSummary/Gene_Centric_Coding_Incl_PTV/"

print(output_path)

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

#######################################################
#     summarize unconditional analysis results
#######################################################

# Combine the lists of results from all 22 chromosomes into a single list named 'results_coding_genome', resulting in a total of 18,435 genes.
# load merged results
results_coding_genome <- get(load(paste0(input_path,"results_coding_genome_batch.Rdata")))

# 7 masks (plof, plof_ds, missense, disruptive_missense, synonymous, ptv, ptv_ds)
results_plof_genome <- NULL
results_plof_ds_genome <- NULL
results_missense_genome <- NULL
results_disruptive_missense_genome <- NULL
results_synonymous_genome <- NULL
results_ptv_genome <- NULL
results_ptv_ds_genome <- NULL

for(kk in 1:length(results_coding_genome))
{
  print(kk)
  
  for (k in 1:7){
    results <- results_coding_genome[[kk]][k][[1]]
    if(is.null(results)==FALSE)
    {
      ### plof
      if(results[3]=="plof")
      {
        results_plof_genome <- rbind(results_plof_genome,results)
      }
      ### plof_ds
      if(results[3]=="plof_ds")
      {
        results_plof_ds_genome <- rbind(results_plof_ds_genome,results)
      }
      ### missense
      if(results[3]=="missense")
      {
        results_missense_genome <- rbind(results_missense_genome,results)
      }
      ### disruptive_missense
      if(results[3]=="disruptive_missense")
      {
        results_disruptive_missense_genome <- rbind(results_disruptive_missense_genome,results)
      }
      ### synonymous
      if(results[3]=="synonymous")
      {
        results_synonymous_genome <- rbind(results_synonymous_genome,results)
      }
      
      ### ptv
      if(results[3]=="ptv")
      {
        results_ptv_genome <- rbind(results_ptv_genome,results)
      }
      ### plof_ds
      if(results[3]=="ptv_ds")
      {
        results_ptv_ds_genome <- rbind(results_ptv_ds_genome,results)
      }
    }
  }
}

###### cMAC_cutoff
# plof
results_plof_genome <- results_plof_genome[results_plof_genome[,"cMAC"]>cMAC_cutoff,]
# plof + disruptive missense
results_plof_ds_genome <- results_plof_ds_genome[results_plof_ds_genome[,"cMAC"]>cMAC_cutoff,]
# missense
results_missense_genome <- results_missense_genome[results_missense_genome[,"cMAC"]>cMAC_cutoff,]
# disruptive missense
results_disruptive_missense_genome <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"cMAC"]>cMAC_cutoff,]
# synonymous
results_synonymous_genome <- results_synonymous_genome[results_synonymous_genome[,"cMAC"]>cMAC_cutoff,]
# ptv
results_ptv_genome <- results_ptv_genome[results_ptv_genome[,"cMAC"]>cMAC_cutoff,]
# ptv + disruptive missense
results_ptv_ds_genome <- results_ptv_ds_genome[results_ptv_ds_genome[,"cMAC"]>cMAC_cutoff,]

### recalculate missense pvalue
if(cMAC_cutoff > 0)
{
  genes_name_disruptive_missense <- as.vector(unlist(results_disruptive_missense_genome[,1])) 
  genes_name_missense <- as.vector(unlist(results_missense_genome[,1])) 
  
  recal_id <- (1:dim(results_missense_genome)[1])[(!genes_name_missense%in%genes_name_disruptive_missense)]
  
  for(kk in recal_id)
  {
    print(kk)
    results_m <- results_missense_genome[kk,]
    
    if(use_SPA)
    {
      ## disruptive missense cMAC < cut_off, set p-value to 1 
      results_missense_genome[kk,"Burden(1,25)-Disruptive"] <- 1
      results_missense_genome[kk,"Burden(1,1)-Disruptive"] <- 1
      
      apc_num <- (length(results_m)-10)/2
      p_seq <- c(1:apc_num,1:apc_num+(apc_num+1))
      
      ## calculate STAAR-B
      pvalues_sub <- as.numeric(results_m[6:length(results_m)][p_seq])
      
      if(sum(is.na(pvalues_sub))>0)
      {
        if(sum(is.na(pvalues_sub))==length(pvalues_sub))
        {
          results_m["STAAR-B"] <- 1
        }else
        {
          ## not all NAs
          pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
          if(sum(pvalues_sub[pvalues_sub<1])>0)
          {
            ## not all ones
            results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
            
          }else
          {
            results_m["STAAR-B"] <- 1
            
          }
        }
      }else
      {
        if(sum(pvalues_sub[pvalues_sub<1])>0)
        {
          results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
        }else
        {
          results_m["STAAR-B"] <- 1
        }
      }
      
      results_missense_genome[kk,"STAAR-B"] <- results_m["STAAR-B"]
      
      ## calculate STAAR-B(1,25)
      pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num)])
      if(sum(is.na(pvalues_sub))>0)
      {
        if(sum(is.na(pvalues_sub))==length(pvalues_sub))
        {
          results_m["STAAR-B(1,25)"] <- 1
        }else
        {
          ## not all NAs
          pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
          if(sum(pvalues_sub[pvalues_sub<1])>0)
          {
            ## not all ones
            results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
            
          }else
          {
            results_m["STAAR-B(1,25)"] <- 1
            
          }
        }
      }else
      {
        if(sum(pvalues_sub[pvalues_sub<1])>0)
        {
          results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
        }else
        {
          results_m["STAAR-B(1,25)"] <- 1
        }
      }
      
      results_missense_genome[kk,"STAAR-B(1,25)"] <- results_m["STAAR-B(1,25)"]
      
      ## calculate STAAR-B(1,1)
      pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1))])
      if(sum(is.na(pvalues_sub))>0)
      {
        if(sum(is.na(pvalues_sub))==length(pvalues_sub))
        {
          results_m["STAAR-B(1,1)"] <- 1
        }else
        {
          ## not all NAs
          pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
          if(sum(pvalues_sub[pvalues_sub<1])>0)
          {
            ## not all ones
            results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
            
          }else
          {
            results_m["STAAR-B(1,1)"] <- 1
            
          }
        }
      }else
      {
        if(sum(pvalues_sub[pvalues_sub<1])>0)
        {
          results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
        }else
        {
          results_m["STAAR-B(1,1)"] <- 1
        }
      }
      
      results_missense_genome[kk,"STAAR-B(1,1)"] <- results_m["STAAR-B(1,1)"]
    }else
    {
      ## disruptive missense cMAC < cut_off, set p-value to 1 
      results_missense_genome[kk,"Burden(1,25)-Disruptive"] <- 1
      results_missense_genome[kk,"Burden(1,1)-Disruptive"] <- 1
      results_missense_genome[kk,"SKAT(1,25)-Disruptive"] <- 1
      results_missense_genome[kk,"SKAT(1,1)-Disruptive"] <- 1			
      results_missense_genome[kk,"ACAT-V(1,25)-Disruptive"] <- 1
      results_missense_genome[kk,"ACAT-V(1,1)-Disruptive"] <- 1			
      
      apc_num <- (length(results_m)-19)/6
      
      p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1))
      
      results_m["STAAR-O"] <- CCT(as.numeric(results_m[6:length(results_m)][p_seq]))
      results_m["STAAR-S(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num)]))
      results_m["STAAR-S(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1))]))
      results_m["STAAR-B(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+2*(apc_num+1))]))
      results_m["STAAR-B(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+3*(apc_num+1))]))
      results_m["STAAR-A(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+4*(apc_num+1))]))
      results_m["STAAR-A(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+5*(apc_num+1))]))
    }
  }
}

###### whole-genome results
# plof
save(results_plof_genome,file="plof.Rdata")
system(paste0("dx upload plof.Rdata --path ",output_path,"plof.Rdata"))
# plof + disruptive missense
save(results_plof_ds_genome,file="plof_ds.Rdata")
system(paste0("dx upload plof_ds.Rdata --path ",output_path,"plof_ds.Rdata"))
# missense
save(results_missense_genome,file="missense.Rdata")
system(paste0("dx upload missense.Rdata --path ",output_path,"missense.Rdata"))
# disruptive missense
save(results_disruptive_missense_genome,file="disruptive_missense.Rdata")
system(paste0("dx upload disruptive_missense.Rdata --path ",output_path,"disruptive_missense.Rdata"))
# synonymous
save(results_synonymous_genome,file="synonymous.Rdata")
system(paste0("dx upload synonymous.Rdata --path ",output_path,"synonymous.Rdata"))
# ptv
save(results_ptv_genome,file="ptv.Rdata")
system(paste0("dx upload ptv.Rdata --path ",output_path,"ptv.Rdata"))
# ptv + disruptive missense
save(results_ptv_ds_genome,file="ptv_ds.Rdata")
system(paste0("dx upload ptv_ds.Rdata --path ",output_path,"ptv_ds.Rdata"))

# If use_SPA=TRUE, use the STAAR-B method to filter significant results; if use_SPA=FALSE, use the STAAR-O method.
if(use_SPA)
{
  ###### significant results
  # plof
  plof_sig <- results_plof_genome[results_plof_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(plof_sig,file="plof_sig.csv")
  system(paste0("dx upload plof_sig.csv --path ",output_path,"plof_sig.csv"))
  
  # missense
  missense_sig <- results_missense_genome[results_missense_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(missense_sig,file="missense_sig.csv")
  system(paste0("dx upload missense_sig.csv --path ",output_path,"missense_sig.csv"))
  
  # synonymous
  synonymous_sig <- results_synonymous_genome[results_synonymous_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(synonymous_sig,file="synonymous_sig.csv")
  system(paste0("dx upload synonymous_sig.csv --path ",output_path,"synonymous_sig.csv"))
  
  # plof_ds
  plof_ds_sig <- results_plof_ds_genome[results_plof_ds_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(plof_ds_sig,file="plof_ds_sig.csv")
  system(paste0("dx upload plof_ds_sig.csv --path ",output_path,"plof_ds_sig.csv"))
  
  # disruptive_missense
  disruptive_missense_sig <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(disruptive_missense_sig,file="disruptive_missense_sig.csv")
  system(paste0("dx upload disruptive_missense_sig.csv --path ",output_path,"disruptive_missense_sig.csv"))
  
  # ptv
  ptv_sig <- results_ptv_genome[results_ptv_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(ptv_sig,file="ptv_sig.csv")
  system(paste0("dx upload ptv_sig.csv --path ",output_path,"ptv_sig.csv"))
  
  # ptv_ds
  ptv_ds_sig <- results_ptv_ds_genome[results_ptv_ds_genome[,"STAAR-B"]<alpha,,drop=FALSE]
  write.csv(ptv_ds_sig,file="ptv_ds_sig.csv")
  system(paste0("dx upload ptv_ds_sig.csv --path ",output_path,"ptv_ds_sig.csv"))
  
  # coding results
  coding_sig <- rbind(plof_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")],
                      missense_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig <- rbind(coding_sig,synonymous_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig <- rbind(coding_sig,plof_ds_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig <- rbind(coding_sig,disruptive_missense_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig <- rbind(coding_sig,ptv_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig <- rbind(coding_sig,ptv_ds_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  
  write.csv(coding_sig,file="coding_sig.csv")
  system(paste0("dx upload coding_sig.csv --path ",output_path,"coding_sig.csv"))
}else
{
  ###### significant results
  # plof
  plof_sig <- results_plof_genome[results_plof_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(plof_sig,file="plof_sig.csv")
  system(paste0("dx upload plof_sig.csv --path ",output_path,"plof_sig.csv"))
  
  # missense
  missense_sig <- results_missense_genome[results_missense_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(missense_sig,file="missense_sig.csv")
  system(paste0("dx upload missense_sig.csv --path ",output_path,"missense_sig.csv"))
  
  # synonymous
  synonymous_sig <- results_synonymous_genome[results_synonymous_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(synonymous_sig,file="synonymous_sig.csv")
  system(paste0("dx upload synonymous_sig.csv --path ",output_path,"synonymous_sig.csv"))
  
  # plof_ds
  plof_ds_sig <- results_plof_ds_genome[results_plof_ds_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(plof_ds_sig,file="plof_ds_sig.csv")
  system(paste0("dx upload plof_ds_sig.csv --path ",output_path,"plof_ds_sig.csv"))
  
  # disruptive_missense
  disruptive_missense_sig <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(disruptive_missense_sig,file="disruptive_missense_sig.csv")
  system(paste0("dx upload disruptive_missense_sig.csv --path ",output_path,"disruptive_missense_sig.csv"))
  
  # ptv
  ptv_sig <- results_ptv_genome[results_ptv_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(ptv_sig,file="ptv_sig.csv")
  system(paste0("dx upload ptv_sig.csv --path ",output_path,"ptv_sig.csv"))
  
  # ptv_ds
  ptv_ds_sig <- results_ptv_ds_genome[results_ptv_ds_genome[,"STAAR-O"]<alpha,,drop=FALSE]
  write.csv(ptv_ds_sig,file="ptv_ds_sig.csv")
  system(paste0("dx upload ptv_ds_sig.csv --path ",output_path,"ptv_ds_sig.csv"))
  
  # coding results
  coding_sig <- rbind(plof_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")],
                      missense_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig <- rbind(coding_sig,synonymous_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig <- rbind(coding_sig,plof_ds_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig <- rbind(coding_sig,disruptive_missense_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig <- rbind(coding_sig,ptv_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig <- rbind(coding_sig,ptv_ds_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  
  write.csv(coding_sig,file="coding_sig.csv")
  system(paste0("dx upload coding_sig.csv --path ",output_path,"coding_sig.csv"))
}

# More relaxed threshold
alpha_1E5 <- 1E-5
if(use_SPA)
{
  ###### significant results
  # plof
  plof_sig_1E5 <- results_plof_genome[results_plof_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(plof_sig_1E5,file="plof_sig_1E5.csv")
  system(paste0("dx upload plof_sig_1E5.csv --path ",output_path,"plof_sig_1E5.csv"))
  
  # missense
  missense_sig_1E5 <- results_missense_genome[results_missense_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(missense_sig_1E5,file="missense_sig_1E5.csv")
  system(paste0("dx upload missense_sig_1E5.csv --path ",output_path,"missense_sig_1E5.csv"))
  
  # synonymous
  synonymous_sig_1E5 <- results_synonymous_genome[results_synonymous_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(synonymous_sig_1E5,file="synonymous_sig_1E5.csv")
  system(paste0("dx upload synonymous_sig_1E5.csv --path ",output_path,"synonymous_sig_1E5.csv"))
  
  # plof_ds
  plof_ds_sig_1E5 <- results_plof_ds_genome[results_plof_ds_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(plof_ds_sig_1E5,file="plof_ds_sig_1E5.csv")
  system(paste0("dx upload plof_ds_sig_1E5.csv --path ",output_path,"plof_ds_sig_1E5.csv"))
  
  # disruptive_missense
  disruptive_missense_sig_1E5 <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(disruptive_missense_sig_1E5,file="disruptive_missense_sig_1E5.csv")
  system(paste0("dx upload disruptive_missense_sig_1E5.csv --path ",output_path,"disruptive_missense_sig_1E5.csv"))
  
  # ptv
  ptv_sig_1E5 <- results_ptv_genome[results_ptv_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(ptv_sig_1E5,file="ptv_sig_1E5.csv")
  system(paste0("dx upload ptv_sig_1E5.csv --path ",output_path,"ptv_sig_1E5.csv"))
  
  # ptv_ds
  ptv_ds_sig_1E5 <- results_ptv_ds_genome[results_ptv_ds_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
  write.csv(ptv_ds_sig_1E5,file="ptv_ds_sig_1E5.csv")
  system(paste0("dx upload ptv_ds_sig_1E5.csv --path ",output_path,"ptv_ds_sig_1E5.csv"))
  
  # coding results
  coding_sig_1E5 <- rbind(plof_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")],
                          missense_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,synonymous_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,plof_ds_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,disruptive_missense_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,ptv_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,ptv_ds_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
  
  write.csv(coding_sig_1E5,file="coding_sig_1E5.csv")
  system(paste0("dx upload coding_sig_1E5.csv --path ",output_path,"coding_sig_1E5.csv"))
}else
{
  ###### significant results
  # plof
  plof_sig_1E5 <- results_plof_genome[results_plof_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(plof_sig_1E5,file="plof_sig_1E5.csv")
  system(paste0("dx upload plof_sig_1E5.csv --path ",output_path,"plof_sig_1E5.csv"))
  
  # missense
  missense_sig_1E5 <- results_missense_genome[results_missense_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(missense_sig_1E5,file="missense_sig_1E5.csv")
  system(paste0("dx upload missense_sig_1E5.csv --path ",output_path,"missense_sig_1E5.csv"))
  
  # synonymous
  synonymous_sig_1E5 <- results_synonymous_genome[results_synonymous_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(synonymous_sig_1E5,file="synonymous_sig_1E5.csv")
  system(paste0("dx upload synonymous_sig_1E5.csv --path ",output_path,"synonymous_sig_1E5.csv"))
  
  # plof_ds
  plof_ds_sig_1E5 <- results_plof_ds_genome[results_plof_ds_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(plof_ds_sig_1E5,file="plof_ds_sig_1E5.csv")
  system(paste0("dx upload plof_ds_sig_1E5.csv --path ",output_path,"plof_ds_sig_1E5.csv"))
  
  # disruptive_missense
  disruptive_missense_sig_1E5 <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(disruptive_missense_sig_1E5,file="disruptive_missense_sig_1E5.csv")
  system(paste0("dx upload disruptive_missense_sig_1E5.csv --path ",output_path,"disruptive_missense_sig_1E5.csv"))
  
  # ptv
  ptv_sig_1E5 <- results_ptv_genome[results_ptv_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(ptv_sig_1E5,file="ptv_sig_1E5.csv")
  system(paste0("dx upload ptv_sig_1E5.csv --path ",output_path,"ptv_sig_1E5.csv"))
  
  # ptv_ds
  ptv_ds_sig_1E5 <- results_ptv_ds_genome[results_ptv_ds_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
  write.csv(ptv_ds_sig_1E5,file="ptv_ds_sig_1E5.csv")
  system(paste0("dx upload ptv_ds_sig_1E5.csv --path ",output_path,"ptv_ds_sig_1E5.csv"))
  
  # coding results
  coding_sig_1E5 <- rbind(plof_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")],
                          missense_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,synonymous_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,plof_ds_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,disruptive_missense_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,ptv_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig_1E5 <- rbind(coding_sig_1E5,ptv_ds_sig_1E5[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  
  write.csv(coding_sig_1E5,file="coding_sig_1E5.csv")
  system(paste0("dx upload coding_sig_1E5.csv --path ",output_path,"coding_sig_1E5.csv"))
}


###########################################################################
#       Manhattan Plot
###########################################################################

# If use_SPA=TRUE, columns are Gene name, Chr, STAAR-B; if use_SPA=FALSE, columns are Gene name, Chr, STAAR-O.

### plof
results_STAAR <- results_plof_genome[,c(1,2,dim(results_plof_genome)[2])]

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

# genes_info
# hgnc_symbol chromosome_name start_position end_position
#      OR2T10               1      248590487    248597700
#       OR2M4               1      248231417    248244679
#     RABGGTB               1       75786197     75795079

genes_info_manhattan[is.na(genes_info_manhattan)] <- 1 # Replace missing values with 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof" # Rename the column "STAAR-B" to "plof"

### plof_ds
results_STAAR <- results_plof_ds_genome[,c(1,2,dim(results_plof_ds_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

# Add plof_ds results
genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof_ds"

### missense
if(!use_SPA)
{
  # If use_SPA=FALSE, columns: Gene name, Chr, STAAR-O
  results_STAAR <- results_missense_genome[,c(1,2,dim(results_missense_genome)[2]-6)]
}else
{
  # If use_SPA=TRUE, columns: Gene name, Chr, STAAR-B
  results_STAAR <- results_missense_genome[,c(1,2,dim(results_missense_genome)[2]-2)]
}

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

# Add missense results
genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "missense"

### disruptive_missense
results_STAAR <- results_disruptive_missense_genome[,c(1,2,dim(results_disruptive_missense_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

# Add disruptive_missense results
genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "disruptive_missense"

### synonymous
results_STAAR <- results_synonymous_genome[,c(1,2,dim(results_synonymous_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

# Add synonymous results
genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "synonymous"

### ptv
results_STAAR <- results_ptv_genome[,c(1,2,dim(results_ptv_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

# Add ptv results
genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "ptv"

### ptv_ds
results_STAAR <- results_ptv_ds_genome[,c(1,2,dim(results_ptv_ds_genome)[2])]

results_m <- c()
for(i in 1:dim(results_STAAR)[2])
{
  results_m <- cbind(results_m,unlist(results_STAAR[,i]))
}

colnames(results_m) <- colnames(results_STAAR)
results_m <- data.frame(results_m,stringsAsFactors = FALSE)
results_m[,2] <- as.numeric(results_m[,2])
results_m[,3] <- as.numeric(results_m[,3])

# Add ptv_ds results
genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "ptv_ds"


## ylim
# The minimum p-value among the gene-centric coding results of the 7 masks (plof, plof_ds, missense, disruptive_missense, synonymous, ptv and ptv_ds).
coding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]])
min_y <- ceiling(-log10(coding_minp)) + 1 #ylim=c(0,min_y)

if(is.infinite(min_y)){
  genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]] <- apply(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]],2,function(x){
    x <- ifelse(x==0,1e-308,x)
    return(x)
  })
  min_y <- 308
}

pch <- c(0,1,2,3,4,5,6)
figure1 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$plof,sig.level=alpha,pch=pch[1],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D"))))

figure2 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$plof_ds,sig.level=alpha,pch=pch[2],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D"))))

figure3 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$missense,sig.level=alpha,pch=pch[3],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D"))))

figure4 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$disruptive_missense,sig.level=alpha,pch=pch[4],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D"))))

figure5 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$synonymous,sig.level=alpha,pch=pch[5],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D"))))

figure6 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$ptv,sig.level=alpha,pch=pch[6],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D"))))

figure7 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$ptv_ds,sig.level=alpha,pch=pch[7],col = c("blue4", "orange3"),ylim=c(0,min_y),
                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D"))))

print("Manhattan plot")

png("gene_centric_coding_manhattan.png", width = 9, height = 6, units = 'in', res = 600)

print(figure1)
print(figure2,newpage = FALSE)
print(figure3,newpage = FALSE)
print(figure4,newpage = FALSE)
print(figure5,newpage = FALSE)
print(figure6,newpage = FALSE)
print(figure7,newpage = FALSE)

dev.off()

system(paste0("dx upload gene_centric_coding_manhattan.png --path ",output_path,"gene_centric_coding_manhattan.png"))

##########################################################
#          Q-Q Plot
##########################################################

print("Q-Q plot")
cex_point <- 1

png("gene_centric_coding_qqplot.png", width = 8, height = 8, units = 'in', res = 600)

### plof
## remove unconverged p-values
observed <- sort(genes_info_manhattan$plof[genes_info_manhattan$plof < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### plof_ds
## remove unconverged p-values
observed <- sort(genes_info_manhattan$plof_ds[genes_info_manhattan$plof_ds < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### missense
## remove unconverged p-values
observed <- sort(genes_info_manhattan$missense[genes_info_manhattan$missense < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### disruptive_missense
## remove unconverged p-values
observed <- sort(genes_info_manhattan$disruptive_missense[genes_info_manhattan$disruptive_missense < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### synonymous
## remove unconverged p-values
observed <- sort(genes_info_manhattan$synonymous[genes_info_manhattan$synonymous < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### ptv
## remove unconverged p-values
observed <- sort(genes_info_manhattan$ptv[genes_info_manhattan$ptv < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=5, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

### ptv_ds
## remove unconverged p-values
observed <- sort(genes_info_manhattan$ptv_ds[genes_info_manhattan$ptv_ds < 1])

lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

par(new=T)
par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=6, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

legend("topleft",legend=c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D"),ncol=1,bty="o",box.lwd=1,pch=0:6,cex=1.5,text.font=2)

dev.off()
system(paste0("dx upload gene_centric_coding_qqplot.png --path ",output_path,"gene_centric_coding_qqplot.png"))

rm(list=setdiff(ls(), c("input_path", "output_path", "alpha", "cMAC_cutoff", "col_c", "use_SPA"))); gc()
print(ls())
