library(tidyverse)
library(scales)
library(showtext)
library(sysfonts)

base_family <- "Times New Roman"

font_add(
  family = "Times New Roman",
  regular = "/usr/share/fonts/TTF/Times.TTF",
  bold = "/usr/share/fonts/TTF/Timesbd.TTF",
  italic = "/usr/share/fonts/TTF/Timesi.TTF",
  bolditalic = "/usr/share/fonts/TTF/Timesbi.TTF"
)

showtext_auto()

project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)
input_dir <- file.path(project_root, "inputs")

infile <- file.path(input_dir, "all_pass_snps_presence_absence_matrix.tsv")
outfile_pdf <- file.path(project_root, "high_low_recurrence.pdf")

mat <- read_tsv(infile, show_col_types = FALSE)

fixed_cols <- c("CHROM", "POS", "REF", "ALT", "variant_id")
sample_cols <- setdiff(colnames(mat), fixed_cols)

low_samples <- c(
  "WRS_1_L",
  "WRS_2_L",
  "WRS_3_L"
)

high_samples <- c(
  "WRS_1_H",
  "WRS_2_H",
  "WRS_3_H"
)

missing_low <- setdiff(low_samples, colnames(mat))
missing_high <- setdiff(high_samples, colnames(mat))

comp_base <- function(x) {
  recode(
    x,
    "A" = "T",
    "T" = "A",
    "C" = "G",
    "G" = "C",
    .default = NA_character_
  )
}

get_mutation_type <- function(ref, alt) {
  ref <- toupper(ref)
  alt <- toupper(alt)

  ref2 <- if_else(ref %in% c("A", "G"), comp_base(ref), ref)
  alt2 <- if_else(ref %in% c("A", "G"), comp_base(alt), alt)

  paste0(ref2, ">", alt2)
}

mut_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

variant_df <- mat %>%
  mutate(
    REF = toupper(REF),
    ALT = toupper(ALT)
  ) %>%
  filter(
    REF %in% c("A", "C", "G", "T"),
    ALT %in% c("A", "C", "G", "T"),
    REF != ALT
  ) %>%
  mutate(
    mutation_type = get_mutation_type(REF, ALT),
    low_degree_raw = rowSums(across(all_of(low_samples)), na.rm = TRUE),
    high_degree_raw = rowSums(across(all_of(high_samples)), na.rm = TRUE),
    low_degree = pmin(low_degree_raw, 3),
    high_degree = pmin(high_degree_raw, 3)
  ) %>%
  filter(mutation_type %in% mut_levels) %>%
  mutate(
    mutation_type = factor(mutation_type, levels = mut_levels)
  )

low_marginal <- variant_df %>%
  filter(low_degree %in% 1:3) %>%
  count(mutation_type, degree = low_degree, name = "count") %>%
  complete(
    mutation_type = mut_levels,
    degree = 1:3,
    fill = list(count = 0)
  ) %>%
  mutate(group = "Low")

high_marginal <- variant_df %>%
  filter(high_degree %in% 1:3) %>%
  count(mutation_type, degree = high_degree, name = "count") %>%
  complete(
    mutation_type = mut_levels,
    degree = 1:3,
    fill = list(count = 0)
  ) %>%
  mutate(group = "High")

marginal_degree <- bind_rows(low_marginal, high_marginal) %>%
  group_by(mutation_type, group) %>%
  mutate(
    total_present_in_group = sum(count),
    fraction = if_else(total_present_in_group > 0, count / total_present_in_group, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    group = factor(group, levels = c("High", "Low")),
    group_short = recode(as.character(group), "High" = "H", "Low" = "L"),
    group_short = factor(group_short, levels = c("H", "L"))
  )

run_degree_test <- function(df_one) {
  tab <- df_one %>%
    select(group, degree, count) %>%
    complete(
      group = c("High", "Low"),
      degree = 1:3,
      fill = list(count = 0)
    ) %>%
    pivot_wider(
      names_from = degree,
      values_from = count
    ) %>%
    arrange(factor(group, levels = c("High", "Low")))

  mat_test <- as.matrix(tab[, c("1", "2", "3")])
  rownames(mat_test) <- tab$group

  if (sum(mat_test) == 0 || any(rowSums(mat_test) == 0)) {
    return(tibble(
      p_value = NA_real_,
      method = NA_character_
    ))
  }

  chi <- suppressWarnings(chisq.test(mat_test))

  tibble(
    p_value = chi$p.value,
    method = "Chi-square test"
  )
}

test_df <- marginal_degree %>%
  group_by(mutation_type) %>%
  group_modify(~ run_degree_test(.x)) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(p_value, method = "holm")
  )

p <- marginal_degree %>%
  ggplot(aes(x = factor(degree), y = fraction, fill = group_short)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "black",
    linewidth = 0.25
  ) +
  facet_wrap(~ mutation_type, nrow = 2) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Number of samples carrying the variant",
    y = "Fraction of variants present in group",
    fill = "Group"
  ) +
  theme_bw(base_size = 14, base_family = base_family) +
  theme(
    text = element_text(family = base_family, colour = "black"),
    strip.text = element_text(size = 14, face = "bold", colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    legend.title = element_text(size = 13, colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(
  outfile_pdf,
  p,
  width = 9,
  height = 5.5,
  device = cairo_pdf
)
