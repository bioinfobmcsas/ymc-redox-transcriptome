library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(showtext)
library(sysfonts)

project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)
input_dir <- file.path(project_root, "inputs")

presence_file <- file.path(input_dir, "all_pass_snps_presence_absence_matrix.tsv")
callable_file <- file.path(input_dir, "callable_base_counts.tsv")
outfile_pdf <- file.path(project_root, "shared_site_subtraction.pdf")

mut_levels <- c(
  "C>T/G>A",
  "C>A/G>T",
  "C>G/G>C",
  "T>A/A>T",
  "T>C/A>G",
  "T>G/A>C"
)

mut_plot_labels <- c(
  "C>T/G>A" = "C>T",
  "C>A/G>T" = "C>A",
  "C>G/G>C" = "C>G",
  "T>A/A>T" = "T>A",
  "T>C/A>G" = "T>C",
  "T>G/A>C" = "T>G"
)

cond_levels <- c("H", "L")
base_family <- "Times New Roman"

font_add(
  family = base_family,
  regular = "/usr/share/fonts/TTF/Times.TTF",
  bold = "/usr/share/fonts/TTF/Timesbd.TTF",
  italic = "/usr/share/fonts/TTF/Timesi.TTF",
  bolditalic = "/usr/share/fonts/TTF/Timesbi.TTF"
)

showtext_auto()

collapse_mut <- function(x) {
  case_when(
    x %in% c("C>A", "G>T") ~ "C>A/G>T",
    x %in% c("C>G", "G>C") ~ "C>G/G>C",
    x %in% c("C>T", "G>A") ~ "C>T/G>A",
    x %in% c("T>A", "A>T") ~ "T>A/A>T",
    x %in% c("T>C", "A>G") ~ "T>C/A>G",
    x %in% c("T>G", "A>C") ~ "T>G/A>C",
    TRUE ~ NA_character_
  )
}

presence <- read.delim(
  presence_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

callable <- read.delim(
  callable_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

sample_cols <- grep("^WRS_[0-9]+_[HL]$", names(presence), value = TRUE)

if (length(sample_cols) == 0) {
  stop("No sample columns matching WRS_<rep>_<H/L> were found.")
}

sample_meta <- tibble(sample = sample_cols) %>%
  mutate(
    rep_id = str_match(sample, "^WRS_([0-9]+)_([HL])$")[, 2],
    cond = str_match(sample, "^WRS_([0-9]+)_([HL])$")[, 3],
    rep_id = factor(rep_id),
    cond = factor(cond, levels = cond_levels)
  )

missing_callable <- setdiff(sample_meta$sample, callable$sample)

if (length(missing_callable) > 0) {
  stop("Samples missing from callable file: ", paste(missing_callable, collapse = ", "))
}

variant_df <- presence %>%
  mutate(
    REF = toupper(REF),
    ALT = toupper(ALT),
    raw_mut = paste0(REF, ">", ALT),
    mut_type = collapse_mut(raw_mut),
    any_H = rowSums(across(matches("^WRS_[0-9]+_H$")), na.rm = TRUE) > 0,
    any_L = rowSums(across(matches("^WRS_[0-9]+_L$")), na.rm = TRUE) > 0,
    phase_shared = any_H & any_L
  ) %>%
  filter(
    REF %in% c("A", "C", "G", "T"),
    ALT %in% c("A", "C", "G", "T"),
    REF != ALT,
    mut_type %in% mut_levels
  ) %>%
  mutate(mut_type = factor(mut_type, levels = mut_levels))

make_rates <- function(df, subtraction_label) {
  counts <- df %>%
    select(mut_type, all_of(sample_cols)) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "present"
    ) %>%
    filter(present > 0) %>%
    count(sample, mut_type, name = "mut_n")

  sample_meta %>%
    expand_grid(mut_type = factor(mut_levels, levels = mut_levels)) %>%
    left_join(counts, by = c("sample", "mut_type")) %>%
    mutate(mut_n = ifelse(is.na(mut_n), 0L, mut_n)) %>%
    left_join(callable, by = "sample") %>%
    mutate(
      A = as.numeric(A),
      C = as.numeric(C),
      G = as.numeric(G),
      T = as.numeric(T),
      callable_n = case_when(
        as.character(mut_type) %in% c("C>A/G>T", "C>G/G>C", "C>T/G>A") ~ C + G,
        as.character(mut_type) %in% c("T>A/A>T", "T>C/A>G", "T>G/A>C") ~ T + A,
        TRUE ~ NA_real_
      ),
      candidate_rate = mut_n / callable_n,
      subtraction = subtraction_label
    ) %>%
    filter(!is.na(callable_n), callable_n > 0)
}

plot_df <- bind_rows(
  make_rates(variant_df, "Before subtraction"),
  make_rates(filter(variant_df, !phase_shared), "After shared-site subtraction")
) %>%
  mutate(
    subtraction = factor(
      subtraction,
      levels = c("Before subtraction", "After shared-site subtraction")
    ),
    mut_type = factor(mut_type, levels = mut_levels)
  )

p <- ggplot(plot_df, aes(x = cond, y = candidate_rate, group = rep_id)) +
  geom_line(color = "grey45", linewidth = 0.45) +
  geom_point(aes(fill = cond), shape = 21, size = 2.8, color = "black", stroke = 0.25) +
  facet_grid(subtraction ~ mut_type, scales = "free_y", labeller = labeller(mut_type = mut_plot_labels)) +
  scale_fill_manual(values = c("H" = "#4DBBD5", "L" = "#E64B35")) +
  labs(
    x = NULL,
    y = "Candidate-site rate",
    fill = "Condition"
  ) +
  theme_bw(base_size = 13, base_family = base_family) +
  theme(
    text = element_text(family = base_family, colour = "black"),
    strip.text = element_text(face = "bold", colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  outfile_pdf,
  p,
  width = 10,
  height = 5.8,
  device = cairo_pdf
)

cat("Done. Results saved to:", outfile_pdf, "\n")
cat("Shared phase candidate sites removed:", sum(variant_df$phase_shared), "\n")
