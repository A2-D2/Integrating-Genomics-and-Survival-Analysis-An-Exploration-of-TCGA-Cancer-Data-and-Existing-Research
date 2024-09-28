#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks", "SummarizedExperiment", "EDASeq")

# Libraries

library(TCGAbiolinks)
library(ggplot2)
library(SummarizedExperiment)

# List the projects in the TCGA Database
 retrieve_tcgaprojects <-getGDCprojects()
 getProjectSummary('TCGA-LUAD')
 

# Retrieve the Lung Adenocarcinoma Data from TCGA
  
 retrieve_luad_data <- GDCquery(project = "TCGA-LUAD",
                                data.category = "Transcriptome Profiling",
                                data.type = "Gene Expression Quantification",
                                experimental.strategy = "RNA-Seq",
                                sample.type = c("Primary Tumor","Solid Tissue Normal"))
 

 # Download the LUAD Data 
 
 GDCdownload(retrieve_luad_data)
 
 # Preparing the LUAD data into an R object
 
 luad_data <- GDCprepare(query=retrieve_luad_data,save=TRUE,save.filename = "luad_data.rda")
 
 # Differentiate between the samples 
 
 # Primary Tumor
 primary_tumor_data <- TCGAquery_SampleTypes(luad_data$barcode,"TP")
 
 # Solid Tissue Normal
 solid_tissue_normal_data <- TCGAquery_SampleTypes(luad_data$barcode,"NT")
 
 
 #Matched with both Primary Tumor and Normal 
 tumor_matched_data <- TCGAquery_MatchedCoupledSampleTypes(c(primary_tumor_data, solid_tissue_normal_data), c("TP","NT"))
 
 # Retrieving the paired data 
 
 retrieve_paired_luad_data <- GDCquery(project = "TCGA-LUAD",
                                data.category = "Transcriptome Profiling",
                                data.type = "Gene Expression Quantification",
                                experimental.strategy = "RNA-Seq",
                                sample.type = c("Primary Tumor","Solid Tissue Normal"),
                                barcode=tumor_matched_data)
 
 GDCdownload(retrieve_paired_luad_data)
 
 
 #Prepare the paired data into an R object
# luad_paired_data <- GDCprepare(retrieve_paired_luad_data,summarizedExperiment = TRUE)
 #luad_paired_data_test <- GDCprepare(query = retrieve_paired_luad_data, )
 luad_paired_data <- GDCprepare(query=retrieve_paired_luad_data,save=TRUE,save.filename = "luad_paired_data.rda")
 
 
 
 
 