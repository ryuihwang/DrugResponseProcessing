##GDSC
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3610/
# Download the dada, ~5Gb
library(readr)
mtx <- read.csv("/Users/hri/gdsc/expression_matrix.csv", row.names = 1, check.names = FALSE)
colnames(mtx) <- sub(pattern = ".cel", "", colnames(mtx), fixed = TRUE)
rownames(mtx) <- sub("^X", "", rownames(mtx))
mtxe <- mtx
sample_annotations <- read_tsv("E-MTAB-3610.sdrf.txt")
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` != '5500994157493061613625_A01', ]
sample_annotations$`Assay Name`
class(colnames(mtx))
colnames(mtx)
common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)
mtxe <- mtx[, colnames(mtx) %in% common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]
all.equal(colnames(mtx), sample_annotations$`Assay Name`)

# Gene annotations
# BiocManager::install("hgu219.db", update = FALSE)
library(hgu219.db)
k <- keys(hgu219.db,keytype="PROBEID")
gene_annotations <- select(hgu219.db, keys=k, columns=c("SYMBOL","GENENAME", "ENSEMBL"), keytype="PROBEID")
gene_annotations <- as_tibble(gene_annotations)
# Get common gene names
common_genes <- intersect(rownames(mtx), gene_annotations$PROBEID)
# Subset and match both matrices
mtx <- mtx[rownames(mtx) %in% common_genes, ]
gene_annotations <- gene_annotations[gene_annotations$PROBEID %in% common_genes, ]
gene_annotations <- gene_annotations[match(rownames(mtx), gene_annotations$PROBEID), ]
all.equal(rownames(mtx), gene_annotations$PROBEID)  # Check if matching worked

library(WGCNA)
is_na_df <- gene_annotations[!complete.cases(gene_annotations), ]
mtx_collapsed <- collapseRows(datET = mtx, rowGroup = gene_annotations$SYMBOL, rowID = rownames(mtx))$datETcollapsed
colnames(mtx_collapsed) <- sample_annotations$`Characteristics[cell line]`
write.csv(mtx_collapsed, file = "final_expression.csv", row.names = TRUE)
