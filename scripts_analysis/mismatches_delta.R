

library(data.table)
library(ggplot2)
library(showtext)
library(sysfonts)
library(ragg)


font_add(

  family = "Times New Roman",
  regular = "/usr/share/fonts/TTF/Times.TTF"

)


showtext_auto()
cond_ref  <- "H"
cond_test <- "L"
project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)
input_dir <- file.path(project_root, "inputs")
res <- fread(file.path(input_dir, "loci_delta.csv"))

scientific_labels <- function(x) {
  labels <- vapply(x, function(value) {
    if (is.na(value)) {
      return("NA")
    }
    if (value == 0) {
      return("0")
    }

    exponent <- floor(log10(abs(value)))
    mantissa <- abs(value) / 10^exponent
    mantissa_label <- format(
      signif(mantissa, 10),
      scientific = FALSE,
      trim = TRUE
    )
    sign_label <- if (value < 0) "−" else ""

    paste0(
      "'", sign_label, mantissa_label, "'",
      " %*% 10^{", exponent, "}"
    )
  }, character(1))

  parse(text = labels)
}


if ("mut_group" %in% names(res) && !"mut_type" %in% names(res)) {
  setnames(res, "mut_group", "mut_type")
}

res <- res[order(-mean_delta)]
res[, mut_type := factor(mut_type, levels = mut_type)]


p_mean <- ggplot(res, aes(x = mut_type, y = mean_delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_col(width = 0.65) +
  geom_errorbar(
    aes(ymin = t_ci_low, ymax = t_ci_high),
    width = 0.15,
    linewidth = 0.8
  ) +
  theme_classic(base_size = 16, base_family = "Times New Roman") +
  scale_y_continuous(labels = scientific_labels) +
  theme(
    text = element_text(family = "Times New Roman", colour = "black"),
    axis.text = element_text(family = "Times New Roman", colour = "black"),
    axis.title = element_text(family = "Times New Roman", colour = "black")
  ) +
  labs(
    x = NULL,
    y = expression("Average " * Delta * "MAF per site")
  )
ggsave(
  file.path(project_root, "mean_delta_by_6_muttypes_ordered.pdf"),
  p_mean,
  width = 7,
  height = 4.8,
  device = cairo_pdf
)
