# DrugResponseProcessing

Overview
This repository contains scripts for downloading, processing, and preparing RNA-Seq and pathology data from various sources, including The Cancer Genome Atlas (TCGA) and other public databases. The processed data is intended for downstream analysis and machine learning applications.

Prerequisites
Ensure the following libraries and packages are installed in your R environment: 
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "hgu219.db", "WGCNA"))
install.packages(c("httr", "data.table", "GenePatternFileReader", "ArrayExpress", "limma", "edgeR", "glmnet", "factoextra", "FactoMineR", "caret", "gplots", "survival", "survminer", "RColorBrewer", "gProfileR", "genefilter", "affy", "preprocessCore", "rhdf5", "dplyr", "stargazer", "readr"))

RNA-Seq Data Processing
1. Query and Download RNA-Seq Data from TCGA:
Use TCGAbiolinks to query and download RNA-Seq data for specific projects.
2. Prepare and Annotate Expression Data:
Convert the expression data to a dataframe, annotate gene IDs, and save the processed data to a CSV file.
3. Download Pathology Images:
Use httr to download pathology image files using UUIDs from a manifest file.
4. Query Clinical Data:
Use TCGAbiolinks to query and download clinical data for specific projects.


Data Integration and Additional Processing
1. Process CCLE Data:
Load and process RNA-Seq data from the Cancer Cell Line Encyclopedia (CCLE).
2. Process GDSC Data:
Load and process expression data from the Genomics of Drug Sensitivity in Cancer (GDSC) project.

Gene Annotations
: Retrieve Gene Annotations:
Use hgu219.db to retrieve and process gene annotations.
