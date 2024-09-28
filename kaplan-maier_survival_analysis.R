# Performing survival analysis on BRCA data from the TCGA Database

# Analysis of TP53 on overall survival 

library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

#getting days of survival data from TCGA

#Clinical BRCA data 

brca_clinical_data <- GDCquery_clinic("TCGA-BRCA")
#
any(colnames(brca_clinical_data) %in% c("vital_status", "days_to_last_follow_up","days_to_death"))
#
which(colnames(brca_clinical_data) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))

brca_clinical_data[,c(9,39,45)]


# Evaluating survival variables
table(brca_clinical_data$vital_status)

# days_to_death is the days between diagnosis until the patient's death
# days_to_last_follow_up is the number of days between the initial diagnosis and the last visit.

# censoring patients that are still alive
brca_clinical_data$deceased <- ifelse(brca_clinical_data$vital_status == "Alive", FALSE, TRUE)

# OS is overall survival. Which is days to death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
brca_clinical_data$OS <- ifelse(brca_clinical_data$vital_status == "Alive",
                                brca_clinical_data$days_to_last_follow_up,
                                brca_clinical_data$days_to_death)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death 
#(clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.


# Retrieving gene expression data using GDCquery
retrieve_brca_data = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")


brca_data <- getResults(retrieve_brca_data)
# get 20 primary tissue sample bar codes
#tumor_data <- brca_data$cases[1:100]
# OR
tumor_data <- brca_data[brca_data$sample_type == "Primary Tumor", "cases"][1:100]




# # get gene expression data from 20 primary tumors 
brca_selected_data <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor_data)

# download data
GDCdownload(brca_selected_data)

# get counts
summarized_experiment_brca <- GDCprepare(brca_selected_data, summarizedExperiment = TRUE)
brca_unstranded_counts <- assay(summarized_experiment_brca, "unstranded")
brca_unstranded_counts[1:10,1:10]


# extract gene and sample metadata from the previously made Summarised Experiment object
gene_data_SE <- as.data.frame(rowData(summarized_experiment_brca))
column_data<- as.data.frame(colData(summarized_experiment_brca))


# Variance Stabilizing Transformation transform counts for survival analysis
# Setting up countData object   
deseq_dataset <- DESeqDataSetFromMatrix(countData = brca_unstranded_counts,
                              colData = column_data,
                              design = ~ 1)


# Removing genes with sum total of 10 reads across all samples
retained_genes <- rowSums(counts(deseq_dataset)) >= 10
deseq_dataset <- deseq_dataset[retained_genes,]


# variance stabilizing transformation 
variance_stabilizing_dataset <- vst(deseq_dataset, blind=FALSE)
brca_vst <- assay(variance_stabilizing_dataset)
brca_vst[1:10,1:10]

# Get data for TP53 gene and add gene metadata information to it -------------
# Using Join to combine the gene metadata and the vst data to get the tp53 data.
tp53_brca_data <- brca_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'cases_id', value = 'count', -gene_id) %>% 
  left_join(., gene_data_SE, by = "gene_id") %>% 
  filter(gene_name == "TP53")


# Retrieve the median value to classify patients into high and low expression counts
median_calculation <- median(tp53_brca_data$count)

# denote which cases have higher or lower expression than median count
tp53_brca_data$strata_data <- ifelse(tp53_brca_data$count >= median_calculation, "HIGH", "LOW")

# Merging clinical data with the TP53 data 
tp53_brca_data$cases_id <- gsub('-01.*', '', tp53_brca_data$cases_id)
tp53_brca_data <- merge(tp53_brca_data, brca_clinical_data, by.x = 'cases_id', by.y = 'submitter_id')

# Computing the survival curve

fitting_curve <- survfit(Surv(OS,deceased) ~ strata_data, data =tp53_brca_data)
fitting_curve

ggsurvplot(fitting_curve,
           data=tp53_brca_data,
           pval = T,
           risk.table = T)

# Using survival diff

fitting_curve_2 <- survdiff(Surv(OS, deceased) ~ strata_data, data = tp53_brca_data)
fitting_curve_2
