###########################################################
# Step 0: Select independent variants from a known variants list to be used in conditional analysis
# Summarize single-variant analysis results and perform conditional analysis of unconditionally significant variants by adjusting a list of known variants.

###########################################################

### Load required R packages
system("dx-restore-folder /STAARpipeline/R_packages/rstudio_workbench_ukbrap_trial.zilinli_iu.2025-01-26T02-24-06.tar.gz")

############## load source code
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

############# Extract CHR, POS, REF and ALT from #rs 

## Info of known variants
# Download from GWAS Catalog
known_loci_sig <- get(load("/mnt/project/UKB_500K_WGS_staarpipeline/CardioAge/Corrected_delta/follow_up_analysis/Conditional_analysis/known_loci_cardiovascular_age_sig.Rdata"))
## rs channel name in aGDS
rs_channel <- "annotation/info/FunctionalAnnotation/rsid"

##### Main Function 

Known_Loci_Info_UKB_RAP <- function(known_loci){
  known_loci_info <- c()
  
  for(chr in 1:22)
  {
    print(chr)
    
    if(sum(known_loci$CHROM==chr)>0){
      
      ## load aGDS files
      gds.name <- paste0("ukb.500k.wgs.chr",chr,".pass.annotated.gds")
      print(gds.name)
      gds.path <- paste0("/mnt/project/UKB_500K_WGS_AGDS_uncompressed_newQClabel_HWE/",gds.name)
      genofile <- seqOpen(gds.path)
      
      variant.id <- seqGetData(genofile, "variant.id")
      rs_num <- seqGetData(genofile,rs_channel)
      
      rs_num_in <- rs_num%in%known_loci$rs
      
      if(sum(rs_num_in) > 0)
      {
        variant.id.in <- variant.id[rs_num_in]
        
        rm(rs_num)
        gc()
        
        rm(variant.id)
        gc()
        
        seqSetFilter(genofile,variant.id=variant.id.in)
        
        ### Basic Info of Significant Loci
        position <- as.numeric(seqGetData(genofile, "position"))
        REF <- as.character(seqGetData(genofile, "$ref"))
        ALT <- as.character(seqGetData(genofile, "$alt"))
        
        known_loci_info_chr <- data.frame(CHR=rep(chr,length(position)),POS=position,REF=REF,ALT=ALT)
        known_loci_info <- rbind(known_loci_info,known_loci_info_chr)	
      }
      seqClose(genofile)
      rm(genofile)
    }
  }
  return(known_loci_info)
}

known_loci_sig_info <- Known_Loci_Info_UKB_RAP(known_loci_sig)

## Output Info of GWAS Catalog SNVs
save(known_loci_sig_info,file="known_loci_sig_info.Rdata")
system("dx upload known_loci_sig_info.Rdata --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/Conditional_analysis/known_loci_sig_info.Rdata")

###########################################################
#           User Input
###########################################################

## Null model
system("dx download -f /UKB_500K_WGS_staarpipeline/CardioAge/revised/WGS_CardioAge_obj.STAAR.UKB.20260125.Rdata")
obj_nullmodel <- get(load("WGS_CardioAge_obj.STAAR.UKB.20260125.Rdata"))

## QC_label
QC_label <- "annotation/info/QC_label3"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## method_cond
method_cond <- "optimal"

## maf_cutoff
samplesize <- length(obj_nullmodel$id_include)
maf_cutoff <- 20.5/samplesize

###########################################################
#           Main Function 
###########################################################

LD_pruning_UKBB_RAP <- function(obj_nullmodel,single_variants_list,maf_cutoff,
                                method_cond="optimal",QC_label="annotation/filter",
                                variant_type="variant",geno_missing_imputation="mean"){
  
  variants_list <- single_variants_list[,c("CHR","POS","REF","ALT")]
  
  known_loci <- c()
  for(chr in 1:22)
  {
    print(chr)
    if(sum(variants_list$CHR==chr)>0)
    {
      variants_list_chr <- variants_list[variants_list$CHR==chr,,drop=FALSE]
      
      ## load aGDS files
      gds.name <- paste0("ukb.500k.wgs.chr",chr,".pass.annotated.gds")
      print(gds.name)
      gds.path <- paste0("/mnt/project/UKB_500K_WGS_AGDS_uncompressed_newQClabel_HWE/",gds.name)
      genofile <- seqOpen(gds.path)
      
      known_loci_chr <- LD_pruning(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,variants_list=variants_list_chr,maf_cutoff=maf_cutoff,
                                   method_cond=method_cond,QC_label=QC_label,
                                   variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)
      
      
      known_loci <- rbind(known_loci,known_loci_chr)
      
      seqClose(genofile)
      rm(genofile)
    }
  }
  
  known_loci_info_annotation <- dplyr::left_join(known_loci,single_variants_list,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))
  return(known_loci_info_annotation)
}


known_loci_sig_LD_pruning <- LD_pruning_UKBB_RAP(obj_nullmodel=obj_nullmodel,
                                                 single_variants_list=known_loci_sig_info,
                                                 maf_cutoff=maf_cutoff,
                                                 method_cond=method_cond,QC_label=QC_label,
                                                 variant_type=variant_type,
                                                 geno_missing_imputation=geno_missing_imputation)

save(known_loci_sig_LD_pruning,file="known_loci_sig_LD_pruning.Rdata")
system("dx upload known_loci_sig_LD_pruning.Rdata --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/Conditional_analysis/known_loci_sig_LD_pruning.Rdata")

write.csv(known_loci_sig_LD_pruning,file ='known_loci_sig_LD_pruning.csv',
          row.names = FALSE)
system("dx upload known_loci_sig_LD_pruning.csv --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/Conditional_analysis/known_loci_sig_LD_pruning.csv")

################# Conditional analysis for Single Variant Analysis 

variant_type <- "variant"
individual_results_sig <- get(load("/mnt/project/UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/results_individual_analysis_sig_5E-8_anno.Rdata"))

Conditional_analysis_UKB_RAP <- function(results_sig,known_loci){
  results_sig_cond <- c()
  for(chr in 1:22){
    if(sum(results_sig$CHR==chr)>=1){
      results_sig_chr <- results_sig[results_sig$CHR==chr,]
      
      ## load aGDS files
      gds.name <- paste0("ukb.500k.wgs.chr",chr,".pass.annotated.gds")
      print(gds.name)
      gds.path <- paste0("/mnt/project/UKB_500K_WGS_AGDS_uncompressed_newQClabel_HWE/",gds.name)
      genofile <- seqOpen(gds.path)
      
      results_sig_cond_chr <- Individual_Analysis_cond(chr=chr,individual_results=results_sig_chr,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                       known_loci=known_loci,method_cond=method_cond,QC_label=QC_label,
                                                       variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)
      
      results_sig_cond <- rbind(results_sig_cond,results_sig_cond_chr)
      
      seqClose(genofile)
      rm(genofile)
    }
  }
  return(results_sig_cond)
}

individual_results_sig_cond <- Conditional_analysis_UKB_RAP(individual_results_sig,known_loci_sig_LD_pruning)

save(individual_results_sig_cond,file="individual_results_sig_cond.Rdata")
system("dx upload individual_results_sig_cond.Rdata --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/Conditional_analysis/individual_results_sig_cond.Rdata")

write.csv(individual_results_sig_cond,file ='individual_results_sig_cond.csv',
          row.names = FALSE)
system("dx upload individual_results_sig_cond.csv --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/Conditional_analysis/individual_results_sig_cond.csv")