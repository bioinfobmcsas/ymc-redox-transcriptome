library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(showtext)

project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)
input_dir <- file.path(project_root, "inputs")

counts_file <- file.path(input_dir, "wrs3_s1_counts.txt")
outfile_pdf <- file.path(project_root, "PCA_featureCounts_Low_High.pdf")

font_add("Times New Roman", regular = "/usr/share/fonts/TTF/Times.TTF")
showtext_auto()

counts <- fread(counts_file, comment.char = "#")

count_cols <- grep("\\.bam$", names(counts), value = TRUE)

count_matrix <- as.matrix(counts[, ..count_cols])
rownames(count_matrix) <- counts$Geneid
colnames(count_matrix) <- sub("\\.bam$", "", colnames(count_matrix))

count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

meta <- data.frame(
  sample = colnames(count_matrix),
  stringsAsFactors = FALSE
) %>%
  mutate(
    replicate = sub("_.*$", "", sample),
    wrs = sub("_dupl$", "", sub("^BRR[0-9]+_", "", sample)),
    condition = case_when(
      wrs == "WRS1" ~ "Low",
      wrs == "WRS2" ~ "High"
    ),
    condition = factor(condition, levels = c("Low", "High"))
  )

rownames(meta) <- meta$sample

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = meta,
  design = ~ condition
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

vst_matrix <- assay(vsd)
vst_matrix <- vst_matrix[apply(vst_matrix, 1, var) > 0, ]

pca <- prcomp(t(vst_matrix), scale. = FALSE)
pca_var <- pca$sdev^2 / sum(pca$sdev^2) * 100

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(meta, by = "sample") %>%
  mutate(label = sub("_dupl$", "", sample))

p <- ggplot(pca_df, aes(PC1, PC2, fill = condition, label = label)) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(
    family = "Times New Roman",
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.8,
    point.padding = 0.8,
    segment.color = "grey60",
    segment.size = 0.3
  ) +
  labs(
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
    fill = "Condition"
  ) +
  scale_fill_manual(values = c("Low" = "#4DBBD5", "High" = "#E64B35")) +
  theme_bw(base_size = 14, base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

ggsave(
  outfile_pdf,
  p,
  width = 6,
  height = 5,
  device = cairo_pdf
)
