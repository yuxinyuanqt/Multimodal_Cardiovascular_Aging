###################################### Whole Genome Sequence (WGS) Analysis of Coding and Noncoding Variants Associated with CardioAG ###################################

###################### Pre-step of association analysis using STAARpipeline
(1) WGS data from 490,549 UK Biobank participants were processed using annotated Genomic Data Structure (aGDS) files generated from project VCF files (UK Biobank Field: 23374) via the vcf2agds toolkit.
    See https://github.com/drarwood/vcf2agds_overview and https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00265-4 for additional details.

(2) A sparse Genetic Relatedness Matrix (GRM) was generated using FastSparseGRM.
    See https://github.com/xihaoli/STAARpipeline-Tutorial and https://github.com/rounakdey/FastSparseGRM for additional details.

###################### Association analysis using STAARpipeline
(1) CardioAge_nullmodel.R:         Step 1: Fit STAAR null model.
(2) STAARpipeline_500K_WGS_CardioAge_batch_mode.txt:         Step 2: Perform WGS association analysis in batch mode using the STAARpipeline app.
(3) STAARpipelineSummary folder:     Summarization and visualization of association analysis results using STAARpipelineSummary.
(4) STAARpipelineSummary_Individual_Analysis_Annotation.R:    Functionally annotate a list of variants.
(5) WGS_Conditional_analysis.R: Perform conditional analysis of unconditionally significant variants by adjusting a list of known variants.

All analyses were conducted on the UK Biobank Research Analysis Platform (RAP):
https://ukbiobank.dnanexus.com/

Please also refer to the STAARpipeline tutorial repository for additional details:
https://github.com/xihaoli/STAARpipeline-Tutorial