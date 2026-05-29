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
input_path <- "/mnt/project/UKB_500K_WGS_staarpipeline/CardioAge/revised/ncRNA/"
output_path <- "/UKB_500K_WGS_staarpipeline/CardioAge/revised/STAARpipelineSummary/ncRNA/"

## results summary
cMAC_cutoff <- 20

## ncRNA_pos
ncRNA_pos <- ncRNA_gene
## alpha level for ncRNA genes
alpha_ncRNA <- 2.5E-06

###########################################################
#           Main Function 
###########################################################

## SPA status
# When there is an imbalance in case-control ratios for binary traits, Saddlepoint Approximations (SPA) may be used to obtain accurate p-values. 
use_SPA <- FALSE

results_ncRNA_genome <- get(load(paste0(input_path,"results_ncRNA_genome_batch.Rdata")))

results_ncRNA_genome <- do.call("rbind",results_ncRNA_genome)

###### cMAC_cutoff
results_ncRNA_genome <- results_ncRNA_genome[results_ncRNA_genome[,"cMAC"]>cMAC_cutoff,]

###### whole-genome results
save(results_ncRNA_genome,file="results_ncRNA_genome.Rdata")
system(paste0("dx upload results_ncRNA_genome.Rdata --path ",output_path,"results_ncRNA_genome.Rdata"))

### ncRNA
if(use_SPA)
{
  ncRNA_sig <- results_ncRNA_genome[results_ncRNA_genome[,"STAAR-B"]<alpha_ncRNA,,drop=FALSE]
}else{
  ncRNA_sig <- results_ncRNA_genome[results_ncRNA_genome[,"STAAR-O"]<alpha_ncRNA,,drop=FALSE]
}
write.csv(ncRNA_sig,file="ncRNA_sig.csv")
system(paste0("dx upload ncRNA_sig.csv --path ",output_path,"ncRNA_sig.csv"))

# More relaxed threshold
alpha_1E5 <- 1E-5
### ncRNA
if(use_SPA)
{
  ncRNA_sig <- results_ncRNA_genome[results_ncRNA_genome[,"STAAR-B"]<alpha_1E5,,drop=FALSE]
}else{
  ncRNA_sig <- results_ncRNA_genome[results_ncRNA_genome[,"STAAR-O"]<alpha_1E5,,drop=FALSE]
}
write.csv(ncRNA_sig,file="ncRNA_sig_1E5.csv")
system(paste0("dx upload ncRNA_sig_1E5.csv --path ",output_path,"ncRNA_sig_1E5.csv"))

## manhattan plot for ncRNA genes
if(!is.null(ncRNA_pos))
{
  results_ncRNA_genome_temp <- data.frame(results_ncRNA_genome[,c(1,2,dim(results_ncRNA_genome)[2])],stringsAsFactors = FALSE)
  results_ncRNA_genome_temp[,2] <- as.numeric(results_ncRNA_genome_temp[,2])
  results_ncRNA_genome_temp[,1] <- as.character(results_ncRNA_genome_temp[,1])
  results_ncRNA_genome_temp[,3] <- as.numeric(results_ncRNA_genome_temp[,3])
  
  ncRNA_gene_pos_results <- dplyr::left_join(ncRNA_pos,results_ncRNA_genome_temp,by=c("chr"="Chr","ncRNA"="Gene.name"))
  ncRNA_gene_pos_results[is.na(ncRNA_gene_pos_results[,5]),5] <- 1
  
  print("ncRNA Manhattan plot")
  
  png("gene_centric_ncRNA_manhattan.png", width = 9, height = 6, units = 'in', res = 600)
  
  print(manhattan_plot(as.numeric(ncRNA_gene_pos_results[,1]), (as.numeric(ncRNA_gene_pos_results[,3])+as.numeric(ncRNA_gene_pos_results[,4]))/2, as.numeric(ncRNA_gene_pos_results[,5]), col = c("blue4", "orange3"),sig.level=alpha_ncRNA))
  
  dev.off()
  system(paste0("dx upload gene_centric_ncRNA_manhattan.png --path ",output_path,"gene_centric_ncRNA_manhattan.png"))
}

## Q-Q plot
cex_point <- 1
observed <- sort(ncRNA_gene_pos_results[ncRNA_gene_pos_results[,5] < 1,5])
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

ncRNA_minp <- min(ncRNA_gene_pos_results[,5])
min_ncRNA_y <- ceiling(-log10(ncRNA_minp)) + 1

print("ncRNA Q-Q plot")
png("gene_centric_ncRNA_qqplot.png", width = 8, height = 8, units = 'in', res = 600)

par(mar=c(5,6,4,4))

plot(lexp,lobs,pch=20, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_ncRNA_y),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

dev.off()
system(paste0("dx upload gene_centric_ncRNA_qqplot.png --path ",output_path,"gene_centric_ncRNA_qqplot.png"))

rm(list=setdiff(ls(), c("input_path", "output_path", "alpha_ncRNA", "cMAC_cutoff", "col_c", "use_SPA"))); gc()
print(ls())

