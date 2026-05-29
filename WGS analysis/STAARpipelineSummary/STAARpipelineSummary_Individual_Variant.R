##########################################################
# Summarization and Visualization of Individual 
# Analysis Results using STAARpipelineSummary
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
input_path <- "/mnt/project/UKB_500K_WGS_staarpipeline/CardioAge/revised/Single/"
output_path <- "/UKB_500K_WGS_staarpipeline/CardioAge/revised/STAARpipelineSummary/Single/"

## alpha level
alpha <- 5E-08 #0.05/(1E07)
## results summary
MAC_cutoff <- 20

###########################################################
#           Main Function 
###########################################################

#######################################################
#     summarize unconditional analysis results
#######################################################

results_single_genome <- get(load(paste0(input_path,"results_single_batch.Rdata")))

# Exclude components from the list that have a class of "try-error".
results_single_genome <- results_single_genome[sapply(results_single_genome, function(x) class(x) != "try-error")]
results_single_genome <- do.call("rbind",results_single_genome) # Convert a list to a data.frame (3737936 obs. of 13 variables)
results_single_genome <- results_single_genome[order(results_single_genome$CHR, results_single_genome$POS),] # sort by CHR and POS

MAC_vector <- 2*(results_single_genome$N)*(results_single_genome$MAF)
results_single_genome <- results_single_genome[MAC_vector>=MAC_cutoff,]

# save results
save(results_single_genome,file="results_single_genome.Rdata")
system(paste0("dx upload results_single_genome.Rdata --path ",output_path,"results_single_genome.Rdata"))

# save significant results (p-value less than the specified alpha)
results_single_genome_sig <- results_single_genome[results_single_genome$pvalue<alpha,,drop=FALSE]

save(results_single_genome_sig,file="results_single_genome_sig.Rdata")
system(paste0("dx upload results_single_genome_sig.Rdata --path ",output_path,"results_single_genome_sig.Rdata"))

## manhattan plot (STAARpipelineSummary R package)
png("manhattan.png", width = 9, height = 6, units = 'in', res = 600)

# Extract the p-values of individual analysis results
pvalue <- results_single_genome$pvalue
CHR <- results_single_genome$CHR
POS <- results_single_genome$POS
rm(results_single_genome)
gc()

# Due to the limitations of floating-point numbers in R, p-values smaller than 1e-308 may be assigned a value of 0. 
# To facilitate the creation of a Manhattan plot, the -log10(p-value) for these loci are replaced with 308.
# https://stackoverflow.com/questions/38165221/r-largest-smallest-representable-numbers
if(min(pvalue)==0)
{
  pvalue_log10 <- -log10(pvalue)
  pvalue_log10[!is.finite(pvalue_log10)] <- 308
  
  print(manhattan_plot(CHR, POS, pvalue_log10, use_logp=TRUE, col = c("blue4", "orange3"),sig.level=alpha))
  
  rm(pvalue_log10)
  
}else
{
  print(manhattan_plot(CHR, POS, pvalue, col = c("blue4", "orange3"),sig.level=alpha))
}

gc()

dev.off()

system(paste0("dx upload manhattan.png --path ",output_path,"manhattan.png"))

## Q-Q plot
observed <- sort(pvalue)
lobs <- -(log10(observed))

if(max(lobs)==Inf)
{
  lobs[lobs==Inf] <- 308
}

# Under the null hypothesis, the expected P-values of SNVs follow a uniform distribution from 0 to 1.
expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

gc()

png("qqplot.png", width = 9, height = 9, units = 'in', res = 600)

par(mar=c(5,6,4,4))
plot(lexp,lobs,pch=20, cex=1, xlim = c(0, max(lexp)), ylim = c(0, max(lobs)),
     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

abline(0, 1, col="red",lwd=2)

dev.off()

system(paste0("dx upload qqplot.png --path ",output_path,"qqplot.png"))