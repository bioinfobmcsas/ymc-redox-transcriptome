library(data.table)
library(topGO)

project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)
input_dir <- file.path(project_root, "inputs")

gaf <- fread(
  file.path(input_dir, "sgd.gaf.gz"),
  sep = "\t",
  header = FALSE,
  skip = "!"
)

setnames(gaf, paste0("V", seq_len(ncol(gaf))))

gaf <- gaf[!(grepl("(^|\\|)NOT(\\||$)", V4))]

gene_col <- "V3"
go_col <- "V5"
aspect_col <- "V9"

make_term2gene_from_gaf_topgo <- function(
    gaf_dt,
    ont = c("BP", "MF", "CC"),
    gene_col = "V3",
    go_col = "V5",
    aspect_col = "V9"
) {
  
  ont <- match.arg(ont)
  
  aspect <- switch(
    ont,
    BP = "P",
    MF = "F",
    CC = "C"
  )
  
  x <- gaf_dt[
    get(aspect_col) == aspect,
    .(
      gene = get(gene_col),
      go = get(go_col)
    )
  ]
  
  x <- x[
    !is.na(gene) & gene != "" &
      !is.na(go) & go != ""
  ]
  
  x <- unique(x)
  
  gene2go <- split(x$go, x$gene)
  gene2go <- lapply(gene2go, unique)
  
  all_genes <- names(gene2go)
  
  gene_vec <- integer(length(all_genes))
  names(gene_vec) <- all_genes
  
  if (length(all_genes) > 0) {
    gene_vec[1] <- 1L
  }
  
  allGenesFactor <- factor(gene_vec, levels = c(0, 1))
  
  GOdata <- new(
    "topGOdata",
    ontology = ont,
    allGenes = allGenesFactor,
    annot = annFUN.gene2GO,
    gene2GO = gene2go
  )
  
  terms <- usedGO(GOdata)
  term2genes_list <- genesInTerm(GOdata, terms)
  
  term2gene <- stack(term2genes_list)
  colnames(term2gene) <- c("gene", "term")
  term2gene <- term2gene[, c("term", "gene")]
  term2gene <- unique(term2gene)
  
  list(
    GOdata = GOdata,
    term2gene = term2gene
  )
}

res_bp <- make_term2gene_from_gaf_topgo(
  gaf,
  ont = "BP",
  gene_col = gene_col,
  go_col = go_col,
  aspect_col = aspect_col
)

res_mf <- make_term2gene_from_gaf_topgo(
  gaf,
  ont = "MF",
  gene_col = gene_col,
  go_col = go_col,
  aspect_col = aspect_col
)

res_cc <- make_term2gene_from_gaf_topgo(
  gaf,
  ont = "CC",
  gene_col = gene_col,
  go_col = go_col,
  aspect_col = aspect_col
)

term2gene_bp <- res_bp$term2gene
term2gene_mf <- res_mf$term2gene
term2gene_cc <- res_cc$term2gene

make_nodup_original <- function(
    term2gene_input,
    minGSSize = 3,
    maxGSSize = 500
) {
  
  term2gene_eff <- copy(term2gene_input)
  setDT(term2gene_eff)
  
  term2gene_eff <- term2gene_eff[, .(
    term = as.character(term),
    gene = as.character(gene)
  )]
  
  term2gene_eff <- term2gene_eff[
    !is.na(term) &
      !is.na(gene)
  ]
  
  term2gene_eff <- unique(
    term2gene_eff,
    by = c("term", "gene")
  )
  
  sizes_eff <- term2gene_eff[
    ,
    .(size = uniqueN(gene)),
    by = term
  ]
  
  keep_terms <- sizes_eff[
    size >= minGSSize & size <= maxGSSize,
    term
  ]
  
  term2gene_eff <- term2gene_eff[
    term %in% keep_terms
  ]
  
  sig_dt <- term2gene_eff[
    order(gene),
    .(
      signature = paste(gene, collapse = "\t"),
      size = .N
    ),
    by = term
  ]
  
  sig_dt <- sig_dt[order(size, term)]
  sig_unique <- sig_dt[!duplicated(signature)]
  
  terms_keep <- sig_unique$term
  
  term2gene_nodup <- term2gene_eff[
    term %in% terms_keep
  ]
  
  return(term2gene_nodup)
}

term2gene_bp_nodup <- make_nodup_original(
  term2gene_bp,
  minGSSize = 3,
  maxGSSize = 500
)

term2gene_mf_nodup <- make_nodup_original(
  term2gene_mf,
  minGSSize = 3,
  maxGSSize = 500
)

term2gene_cc_nodup <- make_nodup_original(
  term2gene_cc,
  minGSSize = 3,
  maxGSSize = 500
)

library(data.table)
library(clusterProfiler)
library(BiocParallel)
library(readxl)
library(dplyr)
library(tidyr)
library(GO.db)
library(AnnotationDbi)
library(openxlsx)
library(ggplot2)
library(stringr)
library(tidytext)
library(showtext)
library(sysfonts)
library(ragg)

font_add(
  family = "Times New Roman",
  regular = "/usr/share/fonts/TTF/Times.TTF",
  bold    = "/usr/share/fonts/TTF/Timesbd.TTF",
  italic  = "/usr/share/fonts/TTF/Timesi.TTF",
  bolditalic = "/usr/share/fonts/TTF/Timesbi.TTF"
)
showtext_auto()

BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)


project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)
input_dir <- file.path(project_root, "inputs")

genes_tbl <- read_excel(file.path(input_dir, "deseq2_LvsH_wrs3_s1_genes.xlsx"))
genes_tbl_full <- as.data.table(genes_tbl)

genes_tbl_full$score <- -sign(genes_tbl_full$log2FoldChange) *
  log10(genes_tbl_full$pvalue)

geneList <- genes_tbl_full$score
names(geneList) <- genes_tbl_full$`gene symbol`

geneList <- geneList[
  !is.na(names(geneList)) &
    names(geneList) != "" &
    !is.na(geneList) &
    is.finite(geneList)
]

geneList <- sort(geneList, decreasing = TRUE)

term2gene_list <- list(
  BP = as.data.table(term2gene_bp),
  MF = as.data.table(term2gene_mf),
  CC = as.data.table(term2gene_cc)
)

ancestor_list <- list(
  BP = AnnotationDbi::as.list(GOBPANCESTOR),
  MF = AnnotationDbi::as.list(GOMFANCESTOR),
  CC = AnnotationDbi::as.list(GOCCANCESTOR)
)

plot_titles <- list(
  BP = "GO Biological Process",
  MF = "GO Molecular Function",
  CC = "GO Cellular Component"
)

run_gsea_go <- function(
    ont,
    term2gene,
    geneList,
    minGSSize = 5,
    maxGSSize = 500,
    p_cutoff = 0.05,
    top_n = 15
) {
  
  term2gene <- as.data.table(term2gene)
  
  term2gene <- term2gene[, .(
    term = as.character(term),
    gene = as.character(gene)
  )]
  
  term2gene <- term2gene[
    !is.na(term) & term != "" &
      !is.na(gene) & gene != ""
  ]
  
  term2gene <- unique(term2gene, by = c("term", "gene"))
  
  term2gene_eff <- term2gene[gene %in% names(geneList)]
  
  sizes_eff <- term2gene_eff[, .(N = uniqueN(gene)), by = term]
  
  keep_terms <- sizes_eff[
    N >= minGSSize & N <= maxGSSize,
    term
  ]
  
  term2gene_eff <- term2gene_eff[term %in% keep_terms]
  
  sig_dt <- term2gene_eff[
    order(gene),
    .(
      signature = paste(gene, collapse = "\t"),
      size = uniqueN(gene)
    ),
    by = term
  ]
  
  dup_sigs <- sig_dt[
    duplicated(signature) |
      duplicated(signature, fromLast = TRUE)
  ]
  
  sig_dt <- sig_dt[order(size, term)]
  sig_unique <- sig_dt[!duplicated(signature)]
  
  terms_keep <- sig_unique$term
  
  term2gene_nodup <- term2gene_eff[
    term %in% terms_keep
  ]
  
  set.seed(42)
  
  gsea <- GSEA(
    geneList = geneList,
    TERM2GENE = term2gene_nodup,
    pvalueCutoff = 1,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  )
  
  res <- as.data.table(gsea@result)
  
  if (nrow(res) == 0) {
    warning("No GSEA results for ", ont)
    return(NULL)
  }
  
  res[, Description := vapply(
    ID,
    function(x) {
      tt <- GOTERM[[x]]
      if (is.null(tt)) {
        x
      } else {
        Term(tt)
      }
    },
    character(1)
  )]
  
  res <- res[order(pvalue)]
  
  sig <- copy(res[p.adjust < p_cutoff])
  sig <- sig[!is.na(ID)]
  
  sig_ids <- sig$ID
  anc_list <- ancestor_list[[ont]]
  
  parents_to_remove <- unique(unlist(lapply(sig_ids, function(go_id) {
    ancestors <- anc_list[[go_id]]
    
    if (is.null(ancestors)) {
      return(character(0))
    }
    
    intersect(ancestors, sig_ids)
  })))
  
  sig_pruned <- sig[!ID %in% parents_to_remove]
  
  sig_pruned[, core_n := lengths(strsplit(core_enrichment, "/"))]
  sig_pruned[, Count := core_n]
  sig_pruned[, GeneRatio := core_n / setSize]
  sig_pruned[, Direction := ifelse(NES > 0, "Low DO", "High DO")]
  
  write.xlsx(
    sig_pruned,
    file = paste0("sig_", tolower(ont), "_pruned.xlsx"),
    rowNames = FALSE
  )
  
  plot_df <- sig_pruned[
    order(Direction, p.adjust),
    .SD[1:min(.N, top_n)],
    by = Direction
  ]
  
  plot_df[, Description_short := str_trunc(Description, width = 55)]
  
  plot_df[, Description_reordered := reorder_within(
    Description_short,
    GeneRatio,
    Direction
  )]
  
  p <- ggplot(plot_df, aes(
    x = GeneRatio,
    y = Description_reordered
  )) +
    geom_point(aes(
      size = Count,
      color = p.adjust
    )) +
    facet_grid(
      Direction ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_y_reordered() +
    scale_size_continuous(
      name = "Count",
      range = c(2, 7)
    ) +
    labs(
      title = plot_titles[[ont]],
      x = "GeneRatio",
      y = NULL,
      color = "FDR"
    ) +
    theme_bw(
      base_size = 12,
      base_family = "Times New Roman"
    ) +
    theme(
      axis.text.y = element_text(
        size = 13,
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        colour = "black"
      ),
      axis.title.x = element_text(
        size = 15,
        colour = "black",
        margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 15,
        colour = "black"
      ),
      strip.text.y = element_text(
        size = 13,
        face = "bold",
        colour = "black"
      ),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        colour = "black",
        size = 17,
        margin = margin(b = 8)
      ),
      legend.title = element_text(
        size = 13,
        colour = "black"
      ),
      legend.text = element_text(
        size = 12,
        colour = "black"
      ),
      legend.key.size = unit(0.6, "cm"),
      panel.grid.minor = element_blank()
    )
      
  ggsave(
    filename = paste0("GO_", ont, ".pdf"),
    plot = p,
    width = 11,
    height = 9,
    device = cairo_pdf
  )
  
  return(list(
    gsea = gsea,
    res = res,
    sig_pruned = sig_pruned,
    plot = p
  ))
}


results <- list()

for (ont in c("BP", "MF", "CC")) {
  results[[ont]] <- run_gsea_go(
    ont = ont,
    term2gene = term2gene_list[[ont]],
    geneList = geneList,
    minGSSize = 5,
    maxGSSize = 500,
    p_cutoff = 0.05,
    top_n = 15
  )
}

