if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("data.table", force = TRUE)
install.packages("GenePatternFileReader")

library("ArrayExpress")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("affy")
library(preprocessCore)
library(rhdf5)
library(dplyr)

install.packages("stargazer")
library(stargazer)

##CCLE
statis <- read.csv("/Users/hri/project/shapiro_results.csv")
statis <- statis[, -1]
rownames(statis) <- statis$Gene
statis$Gene <- NULL
stargazer(statis,title="Trained Data Distribution", type="html", out = "stats.html", flip=T)


data <- h5read("/Users/hri/gdsc/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "dataset_name")
AA <- read.delim(file="/Users/hri/gdsc/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz", check.names = FALSE)
AA$gene_id <- sub("\\..*", "", AA$gene_id)
AA <- subset(AA, select = -transcript_ids)

result <- left_join(AA, gene_annotations[, c("ENSEMBL", "SYMBOL")], by = c("gene_id" = "ENSEMBL"))
result <- select(result, ncol(result), 1:(ncol(result)-1))
result <- result[, c(ncol(result), 1:(ncol(result)-1))]
result <- distinct(result, SYMBOL, .keep_all = TRUE)
result <- subset(result, select = -gene_id)
result <- result[complete.cases(result), ]
rownames(result) <- result[, 1]
result <- result[, -1]

for (i in 1:nrow(line)) {
  colnames(result)[colnames(result) == line$CCLE_ID[i]] <- line$Name[i]
}


line <- read.deliresultline <- read.delim(file="/Users/hri/gdsc/Cell_lines_annotations_20181226.txt")
gene <- read_delim(file="/Users/hri/gdsc/gencode.v19.genes.v7_model.patched_contigs.gtf.gz", col_names = TRUE, skip = 6)
gtf_df <- read_delim(gtf_file, delim = "\t", quote = "", comment = "#", col_names = FALSE)

write.csv(result, file = "CCLE_RNA_seq_TPM.csv", row.names = TRUE)