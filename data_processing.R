
# Get the paired data into a single R object for data processing 
paired_summarized_experiment<-get(load("/Users/abhi/lung_adenocarcinoma_analysis/paired_data/luad_paired_data.rda"))

# Now remove the donors with low tumor percentage 
low_tumor_purity<- TCGAtumor_purity(colnames(paired_summarized_experiment),0,0,0,0,0.7) #32
low_tumor_purity_2 <- TCGAtumor_purity(colnames(paired_summarized_experiment),0,0,0,0,0.5) # 47


# Remaining pure tumour barcodes
length(low_tumor_purity$pure_barcodes)

# Donors to be filtered out
length(low_tumor_purity$filtered)

paired_summarized_experiment <- paired_summarized_experiment[,colnames(paired_summarized_experiment)
                                                             %in% union(low_tumor_purity$pure_barcodes,low_tumor_purity$filtered)]
# Now filter the samples from the paired donors

retrieve_paired <- TCGAquery_MatchedCoupledSampleTypes(colnames(paired_summarized_experiment),c("NT","TP"))

paired_summarized_experiment <- paired_summarized_experiment[,colnames(paired_summarized_experiment) %in% retrieve_paired]
save(paired_summarized_experiment,file = "/Users/abhi/lung_adenocarcinoma_analysis/paired_data/luad_paired_highTumorPurity.rda")

# Data Pre-processing
data_preprocessing <-TCGAanalyze_Preprocessing(object = paired_summarized_experiment,cor.cut = 0.7)

# QA - Quality Assurance of the Data as it is still raw data downloaded from the TCGA database
# Removing unexpressed genes 

is_expressed <- data_preprocessing >=10

expression_data_frame <- data.frame(Expressed = rowSums(is_expressed))
ggplot(expression_data_frame,aes(x=Expressed)) + geom_bar()

# Retaining the genes that are expressed in 10 of the samples
retained_genes <- rowSums(data_preprocessing >=10) >=10
table(retained_genes)

# Removing genes that have less than 10 count 

data_preprocessing <- data_preprocessing[retained_genes,]

# Visualizing the distribution of the samples
boxplot(data_preprocessing)

#Making the data easier to view by performing log10
boxplot(log10(data_preprocessing))

normalize_data <- TCGAanalyze_Normalization(tabDF = data_preprocessing,
                                      geneInfo = geneInfoHT, # Using HT rather than usual. explain this
                                      method = "gcContent") 

filter_data <- TCGAanalyze_Filtering(tabDF = normalize_data,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

save(filter_data, file = "/Users/abhi/lung_adenocarcinoma_analysis/paired_data/luad_preprocessed_highTumorPurity.rda")
