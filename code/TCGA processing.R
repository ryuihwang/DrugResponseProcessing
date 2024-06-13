##TCGA_etc
stadquery <- GDCquery(project = "TCGA-OV", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts",
                      experimental.strategy = "RNA-Seq")
GDCdownload(query = stadquery, method = "api",)
stadprpr <- GDCprepare(query = stadquery, summarizedExperiment = T)
assay(stadprpr, "fpkm_uq_unstrand")
colData(stadprpr)
rowData(stadprpr)
expression_HNSC <- as.data.frame(assay(stadprpr, "fpkm_uq_unstrand"))
expression_df <- cbind(expression_df, expression_HNSC)

a <- rownames(expression_df)
expression_df$gene_id <- a
expr <- expression_df[, c(ncol(expression_df), 1:(ncol(expression_df)-1))]
expr <- expr[!grepl("_", expr$gene_id), ]
expr$gene_id <- sub("\\..*", "", expr$gene_id)
expression_df <- left_join(expr, gene_annotations[, c("ENSEMBL", "SYMBOL")], by = c("gene_id" = "ENSEMBL"))
expression_df <- expression_df[, c(ncol(expression_df), 1:(ncol(expression_df)-1))]
expression_df <- distinct(expression_df, SYMBOL, .keep_all = TRUE)
expression_df <- subset(expression_df, select = -gene_id)
expression_df <- expression_df[complete.cases(expression_df), ]
rownames(expression_df) <- expression_df[, 1]
expression_df <- expression_df[, -1]

write.csv(expression_df, file = "TCGA_RNA_seq_FPKM_up.csv", row.names = TRUE)

##pathology data
library(httr)

manifest <- read.table("gdc_manifest.2024-04-07.txt", header = TRUE, stringsAsFactors = FALSE)

download_file <- function(uuid, filename) {
  url <- paste0("https://api.gdc.cancer.gov/files/", uuid)
  response <- httr::GET(url, httr::write_disk(filename, overwrite = TRUE))
  if (httr::http_status(response)$category == "success") {
    cat("Downloaded file:", filename, "\n")
  } else {
    cat("Failed to download file:", filename, "\n")
  }
}

for (i in 1:nrow(manifest)) {
  download_file(manifest$id[i], manifest$md5[i])
}


clin <- GDCquery_clinic("TCGA-LUAD", type = "clinical", save.csv = TRUE)
clin <- GDCquery_clinic("TCGA-ACC", type = "biospecimen", save.csv = TRUE)
clin.cptac2 <- GDCquery_clinic("CPTAC-2", type = "clinical")
clin.TARGET_ALL_P1 <- GDCquery_clinic("TARGET-ALL-P1", type = "clinical")
clin.fm_ad <- GDCquery_clinic("FM-AD", type = "clinical")
## Not run: 
clin <- GDCquery_clinic(project = "CPTAC-3", type = "clinical")
clin <- GDCquery_clinic(project = "CPTAC-2", type = "clinical")
clin <- GDCquery_clinic(project = "HCMI-CMDC", type = "clinical")
clin <- GDCquery_clinic(project = "NCICCR-DLBCL", type = "clinical")
clin <- GDCquery_clinic(project = "ORGANOID-PANCREATIC", type = "clinical")