#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Sc.sgd.db)
})

## -------- CONFIG --------
project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)

COUNTS_FILE <- Sys.getenv(
  "DESEQ2_COUNTS_FILE",
  unset = file.path(project_root, "inputs", "wrs3_s1_counts.txt")
)
SGD_GAF_FILE <- Sys.getenv(
  "DESEQ2_GAF_FILE",
  unset = file.path(project_root, "inputs", "sgd.gaf.gz")
)
OUT_FILE <- Sys.getenv(
  "DESEQ2_OUTPUT_FILE",
  unset = file.path(project_root, "inputs", "deseq2_LvsH_wrs3_s1_genes.xlsx")
)

samples <- c(
  "BRR1_WRS1", "BRR1_WRS2",
  "BRR2_WRS1", "BRR2_WRS2",
  "BRR3_WRS1", "BRR3_WRS2"
)
count_columns <- paste0(samples, "_dupl.bam")

conditions <- factor(
  c("low", "high", "low", "high", "low", "high"),
  levels = c("high", "low")
)
experiments <- factor(c("BRR1", "BRR1", "BRR2", "BRR2", "BRR3", "BRR3"))

## -------- LOAD COUNTS --------
raw <- fread(COUNTS_FILE, sep = "\t", header = TRUE, comment.char = "#")

if (!"Geneid" %in% names(raw)) {
  stop("Column 'Geneid' not found in counts file: ", COUNTS_FILE)
}

missing_samples <- setdiff(count_columns, names(raw))
if (length(missing_samples) > 0) {
  stop("These samples are not in the count matrix: ", paste(missing_samples, collapse = ", "))
}

count_df <- as.data.frame(raw[, ..count_columns])
rownames(count_df) <- raw$Geneid

## -------- RUN DESEQ2 --------
coldata <- data.frame(
  sample = samples,
  condition = conditions,
  experiment = experiments,
  row.names = count_columns
)

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_df)),
  colData = coldata,
  design = ~ experiment + condition
)
dds <- DESeq(dds, quiet = TRUE)

res <- results(
  dds,
  contrast = c("condition", "low", "high"),
  alpha = 0.05
)

res_dt <- as.data.table(as.data.frame(res), keep.rownames = "Geneid")

## -------- BUILD CURRENT GAF ORF -> SYMBOL MAP --------
gaf <- fread(
  SGD_GAF_FILE,
  sep = "\t",
  header = FALSE,
  comment.char = "!",
  fill = TRUE,
  quote = ""
)

if (ncol(gaf) < 12) {
  stop("GAF file has fewer than 12 columns: ", SGD_GAF_FILE)
}


orf_to_sgd <- as.data.table(AnnotationDbi::select(
  org.Sc.sgd.db,
  keys = res_dt$Geneid,
  keytype = "ORF",
  columns = "SGD"
))
setnames(orf_to_sgd, c("ORF", "SGD"), c("Geneid", "SGD_ID"))
orf_to_sgd <- unique(orf_to_sgd[!is.na(SGD_ID) & SGD_ID != ""])

ambiguous_orfs <- orf_to_sgd[, .(n_ids = uniqueN(SGD_ID)), by = Geneid][n_ids > 1]
if (nrow(ambiguous_orfs) > 0) {
  stop(
    "org.Sc.sgd.db maps some ORFs to multiple SGD IDs: ",
    paste(ambiguous_orfs$Geneid, collapse = ", ")
  )
}

sgd_to_symbol <- unique(gaf[V12 != "protein_complex", .(
  SGD_ID = V2,
  `gene symbol` = V3
)])

ambiguous_ids <- sgd_to_symbol[, .(n_symbols = uniqueN(`gene symbol`)), by = SGD_ID][n_symbols > 1]
if (nrow(ambiguous_ids) > 0) {
  stop(
    "Current GAF maps some SGD IDs to multiple symbols: ",
    paste(ambiguous_ids$SGD_ID, collapse = ", ")
  )
}

gaf_map <- merge(orf_to_sgd, sgd_to_symbol, by = "SGD_ID", all.x = TRUE)
gaf_map <- unique(gaf_map[, .(Geneid, `gene symbol`)])

## -------- ADD CURRENT SYMBOLS AND WRITE XLSX --------
res_dt <- merge(res_dt, gaf_map, by = "Geneid", all.x = TRUE, sort = FALSE)
res_dt[is.na(`gene symbol`) | `gene symbol` == "", `gene symbol` := Geneid]
res_dt <- res_dt[order(pvalue, na.last = TRUE)]
setcolorder(res_dt, c(
  "Geneid", "baseMean", "log2FoldChange", "lfcSE",
  "stat", "pvalue", "padj", "gene symbol"
))

wb <- createWorkbook()
addWorksheet(wb, "deseq2_LvsH_wrs3_s1")
writeData(wb, "deseq2_LvsH_wrs3_s1", res_dt)
freezePane(wb, "deseq2_LvsH_wrs3_s1", firstRow = TRUE)
saveWorkbook(wb, OUT_FILE, overwrite = TRUE)

cat("Wrote", nrow(res_dt), "DESeq2 rows with current GAF symbols to", OUT_FILE, "\n")
