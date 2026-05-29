#####################################################
# Functionally annotate a list of variants
# 07/23/2024
#####################################################

### Load required R packages
system("dx-restore-folder /STAARpipeline/R_packages/rstudio_workbench_ukbrap_trial.zilinli_iu.2025-01-26T02-24-06.tar.gz")

############## load source code
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)

################ Input 

## single variant analysis results: a data frame containing the information of variants to be functionally annotated. The data frame must include 4 columns with
#' the following names: "CHR" (chromosome number), "POS" (position), "REF" (reference allele), and "ALT" (alternative allele).
individual_results_sig <- get(load("/mnt/project/UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/results_single_genome_sig_5E-8.Rdata"))

## QC_label: channel name of the QC label in the GDS/aGDS file  (default = "annotation/filter")
QC_label <- "annotation/info/QC_label3"

## Annotation_dir: channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation")
Annotation_dir <- "annotation/info/FunctionalAnnotation"

## Annotation channel
system("dx download -f Annotation_name_catalog.csv")
Annotation_name_catalog <- read.csv("Annotation_name_catalog.csv")

## Use_annotation_weights
Use_annotation_weights <- TRUE
Annotation_name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")


################# main analysis

Annotate_Single_Variants_UKBB_RAP <- function(single_variants_list,
                                              QC_label="annotation/filter",Annotation_dir="annotation/info/FunctionalAnnotation",
                                              Annotation_name_catalog,Annotation_name){
  
  single_variants_list_info <- single_variants_list[,c("CHR","POS","REF","ALT")]
  
  single_variants_list_annotation <- c()
  for(chr in 1:22)
  {
    print(chr)
    if(sum(single_variants_list_info$CHR==chr)>0)
    {
      single_variants_list_info_chr <- single_variants_list_info[single_variants_list_info$CHR==chr,,drop=FALSE]
      
      ## load aGDS files
      gds.name <- paste0("ukb.500k.wgs.chr",chr,".pass.annotated.gds")
      print(gds.name)
      gds.path <- paste0("/mnt/project/UKB_500K_WGS_AGDS_uncompressed_newQClabel_HWE/",gds.name)
      genofile <- seqOpen(gds.path)
      
      position <- as.numeric(seqGetData(genofile, "position"))
      REF <- as.character(seqGetData(genofile, "$ref"))
      ALT <- as.character(seqGetData(genofile, "$alt"))
      variant_id <- seqGetData(genofile, "variant.id")
      
      chr_info <- data.frame(CHR=rep(chr,length(position)),POS=position,REF=REF,ALT=ALT,variant_id=variant_id)
      
      single_variants_list_info_chr <- dplyr::left_join(single_variants_list_info_chr,chr_info,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))
      variant.id.in <- single_variants_list_info_chr$variant_id[!is.na(single_variants_list_info_chr$variant_id)]
      
      seqSetFilter(genofile,variant.id=variant.id.in)
      
      CHR <- as.numeric(seqGetData(genofile, "chromosome"))
      position <- as.numeric(seqGetData(genofile, "position"))
      REF <- as.character(seqGetData(genofile, "$ref"))
      ALT <- as.character(seqGetData(genofile, "$alt"))
      filter <- seqGetData(genofile, QC_label)
      
      Anno.Int.PHRED.sub <- NULL
      Anno.Int.PHRED.sub.name <- NULL
      
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      
      single_variants_list_annotation_chr <- data.frame(CHR=CHR,POS=position,REF=REF,ALT=ALT,QC_label=filter)
      single_variants_list_annotation_chr <- cbind(single_variants_list_annotation_chr,Anno.Int.PHRED.sub)
      
      single_variants_list_annotation <- rbind(single_variants_list_annotation,single_variants_list_annotation_chr)
      seqClose(genofile)
      # system("rm *.gds") # remove aGDS file
    }
    
  }
  single_variants_list_info_annotation <- dplyr::left_join(single_variants_list,single_variants_list_annotation,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))
  
  return(single_variants_list_info_annotation)
}

results_individual_analysis_sig_anno <- Annotate_Single_Variants_UKBB_RAP(single_variants_list=individual_results_sig,
                                                                          QC_label=QC_label,Annotation_dir=Annotation_dir,
                                                                          Annotation_name_catalog=Annotation_name_catalog,
                                                                          Annotation_name=Annotation_name)

save(results_individual_analysis_sig_anno,file="results_individual_analysis_sig_anno.Rdata")
system("dx upload results_individual_analysis_sig_anno.Rdata --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/results_individual_analysis_sig_5E-8_anno.Rdata")

write.csv(results_individual_analysis_sig_anno,file ='results_individual_analysis_sig_anno.csv',
          row.names = FALSE)
system("dx upload results_individual_analysis_sig_anno.csv --path /UKB_500K_WGS_staarpipeline/CardioAge/revised/follow_up_analysis/results_individual_analysis_sig_5E-8_anno.csv")