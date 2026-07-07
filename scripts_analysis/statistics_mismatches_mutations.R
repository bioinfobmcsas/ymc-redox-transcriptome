

library(data.table)

library(dplyr)

library(tidyr)

library(glmmTMB)

library(emmeans)

library(ggplot2)

library(patchwork)

library(showtext)

library(sysfonts)


mut_levels <- c(

  "C>T/G>A",

  "C>A/G>T",

  "C>G/G>C",

  "T>A/A>T",

  "T>C/A>G",

  "T>G/A>C"

)


mut_label <- paste(mut_levels, collapse = " + ")


mut_plot_labels <- c(

  "C>T/G>A" = "C>T",

  "C>A/G>T" = "C>A",

  "C>G/G>C" = "C>G",

  "T>A/A>T" = "T>A",

  "T>C/A>G" = "T>C",

  "T>G/A>C" = "T>G"

)


project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)

input_dir <- file.path(project_root, "inputs")

out_dir <- project_root

collapse_existing_mut_type <- function(x) {

  data.table::fcase(

    x == "C>A", "C>A/G>T",

    x == "C>G", "C>G/G>C",

    x == "C>T", "C>T/G>A",

    x == "T>A", "T>A/A>T",

    x == "T>C", "T>C/A>G",

    x == "T>G", "T>G/A>C",

    default = x

  )

}


infile <- file.path(input_dir, "loci_per_sample_alt_fraction.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


cond_ref  <- "H"

cond_test <- "L"


cond_levels <- c(cond_ref, cond_test)


base_family <- "Times New Roman"


font_add(

  family = "Times New Roman",

  regular = "/usr/share/fonts/TTF/Times.TTF",

  bold    = "/usr/share/fonts/TTF/Timesbd.TTF",

  italic  = "/usr/share/fonts/TTF/Timesi.TTF",

  bolditalic = "/usr/share/fonts/TTF/Timesbi.TTF"

)


showtext_auto()


df <- fread(infile)


df[, mut_type := collapse_existing_mut_type(as.character(mut_type))]


required_cols <- c(

  "sample",

  "replicate",

  "cond",

  "mut_type",

  "COUNT_REF",

  "COUNT_ALT",

  "total",

  "alt_fraction"

)


missing_cols <- setdiff(required_cols, names(df))


if (length(missing_cols) > 0) {

  stop(

    "Missing required columns: ",

    paste(missing_cols, collapse = ", ")

  )

}


for (cc in c("COUNT_REF", "COUNT_ALT", "total", "alt_fraction")) {

  df[, (cc) := as.numeric(get(cc))]

}


df_target <- df %>%

  filter(

    !is.na(sample),

    !is.na(replicate),

    !is.na(cond),

    !is.na(mut_type),

    !is.na(COUNT_REF),

    !is.na(COUNT_ALT),

    COUNT_REF >= 0,

    COUNT_ALT >= 0,

    COUNT_REF + COUNT_ALT > 0

  ) %>%

  mutate(

    mut_type = as.character(mut_type)

  ) %>%

  filter(

    mut_type %in% mut_levels,

    cond %in% cond_levels

  ) %>%

  group_by(sample, replicate, cond, mut = mut_type) %>%

  summarise(

    COUNT_REF = sum(COUNT_REF, na.rm = TRUE),

    COUNT_ALT = sum(COUNT_ALT, na.rm = TRUE),

    .groups = "drop"

  ) %>%

  mutate(

    total = COUNT_REF + COUNT_ALT,

    alt_fraction = COUNT_ALT / total,

    mut = factor(mut, levels = mut_levels),

    cond = factor(cond, levels = cond_levels),

    replicate = factor(replicate)

  ) %>%

  filter(

    total > 0,

    !is.na(mut),

    !is.na(cond),

    !is.na(replicate)

  )


if (nrow(df_target) == 0) {

  stop("No target rows left after filtering: ", mut_label)

}


fit_bb <- glmmTMB(

  cbind(COUNT_ALT, COUNT_REF) ~ cond * mut + (1 | replicate),

  data = df_target,

  family = betabinomial(link = "logit"),

  control = glmmTMBControl(

    optCtrl = list(iter.max = 1e5, eval.max = 1e5)

  )

)


emm_bb <- emmeans(

  fit_bb,

  ~ cond | mut,

  type = "response"

)


emm_bb_dt <- as.data.frame(summary(emm_bb)) %>%

  mutate(

    mut = factor(mut, levels = mut_levels),

    cond = factor(cond, levels = cond_levels)

  ) %>%

  arrange(mut, cond)


pval_to_stars <- function(p) {

  case_when(

    is.na(p) ~ "",

    p < 0.001 ~ "***",

    p < 0.01  ~ "**",

    p < 0.05  ~ "*",

    TRUE ~ ""

  )

}


contrast_bb <- contrast(

  emm_bb,

  method = "revpairwise",

  by = "mut",

  type = "response"

)


posthoc_bb <- as.data.frame(summary(

  contrast_bb,

  infer = TRUE,

  adjust = "none"

)) %>%

  mutate(

    mut = factor(mut, levels = mut_levels)

  ) %>%

  arrange(mut) %>%

  rename(p_value = p.value) %>%

  mutate(

    p_holm = p.adjust(p_value, method = "holm"),

    p_label = pval_to_stars(p_holm)

  ) %>%

  select(

    any_of(c(

      "contrast", "mut", "odds.ratio", "ratio", "SE", "df",

      "asymp.LCL", "asymp.UCL", "null", "z.ratio", "t.ratio",

      "p_value", "p_holm", "p_label"

    )),

    everything()

  )


sig_dt <- posthoc_bb %>%

  select(mut, p_value, p_holm, p_label)


posthoc_fc_dt <- posthoc_bb


if ("odds.ratio" %in% names(posthoc_fc_dt)) {

  posthoc_fc_dt <- posthoc_fc_dt %>%

    rename(MR_fold_change = odds.ratio)

} else if ("ratio" %in% names(posthoc_fc_dt)) {

  posthoc_fc_dt <- posthoc_fc_dt %>%

    rename(MR_fold_change = ratio)

} else {

  stop("Could not find odds.ratio/ratio column in posthoc beta-binomial output.")

}


posthoc_fc_dt <- posthoc_fc_dt %>%

  mutate(mut = factor(mut, levels = mut_levels)) %>%

  arrange(mut)


sig_pos_A <- emm_bb_dt %>%

  group_by(mut) %>%

  summarise(

    y_pos = max(asymp.UCL, na.rm = TRUE) * 1.12,

    .groups = "drop"

  ) %>%

  left_join(sig_dt, by = "mut") %>%

  mutate(mut = factor(mut, levels = mut_levels))


p_est_bb_bar <- ggplot(

  emm_bb_dt,

  aes(x = mut, y = prob, fill = cond)

) +

  geom_col(

    position = position_dodge(width = 0.8),

    width = 0.7,

    color = "black",

    linewidth = 0.2

  ) +

  geom_errorbar(

    aes(ymin = asymp.LCL, ymax = asymp.UCL),

    position = position_dodge(width = 0.8),

    width = 0.2,

    linewidth = 0.4

  ) +

  geom_text(

    data = sig_pos_A,

    aes(x = mut, y = y_pos, label = p_label),

    inherit.aes = FALSE,

    size = 5,

    family = base_family

  ) +

  scale_x_discrete(labels = mut_plot_labels, drop = FALSE) +

  scale_y_continuous(

    expand = expansion(mult = c(0, 0.18))

  ) +

  theme_classic(base_size = 15, base_family = base_family) +

  labs(

    x = "Substitution",

    y = "Estimated MR",

    fill = "Condition"

  ) +

  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))


p_fc_bb <- ggplot(

  posthoc_fc_dt,

  aes(x = mut, y = MR_fold_change)

) +

  geom_hline(yintercept = 1, linetype = "dashed") +

  geom_col(

    width = 0.7,

    color = "black",

    linewidth = 0.2

  ) +

  geom_errorbar(

    aes(ymin = asymp.LCL, ymax = asymp.UCL),

    width = 0.2,

    linewidth = 0.4

  ) +

  geom_text(

    aes(

      y = asymp.UCL * 1.05,

      label = p_label

    ),

    size = 5,

    family = base_family

  ) +

  scale_x_discrete(labels = mut_plot_labels, drop = FALSE) +

  scale_y_continuous(

    limits = c(0, NA),

    expand = expansion(mult = c(0, 0.18))

  ) +

  theme_classic(base_size = 15, base_family = base_family) +

  labs(

    x = "Substitution",

    y = "Estimated MR fold change\nLow DO / High DO"

  ) +

  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))


p_est_bb_bar_clean <- p_est_bb_bar +

  labs(tag = "A") +

  theme(

    legend.position = "right",

    plot.tag = element_text(size = 18, face = "bold", family = base_family),

    plot.tag.position = c(-0.03, 1.03),

    plot.margin = margin(t = 12, r = 8, b = 5, l = 18)

  )


p_fc_bb_clean <- p_fc_bb +

  labs(tag = "B") +

  theme(

    legend.position = "none",

    plot.tag = element_text(size = 18, face = "bold", family = base_family),

    plot.tag.position = c(-0.03, 1.03),

    plot.margin = margin(t = 12, r = 8, b = 5, l = 18)

  )


p_combined_bb <- p_est_bb_bar_clean + p_fc_bb_clean +

  plot_layout(widths = c(1.25, 1))


ggsave(

  file.path(out_dir, "betabinomial_estimated_MR_and_fold_change.pdf"),

  p_combined_bb,

  width = 13,

  height = 5

)


cat("Done.\n")

cat("Results saved to:", out_dir, "\n")

cat("Target substitutions:", mut_label, "\n")

cat("Rows used in model:", nrow(df_target), "\n")

for (m in mut_levels) {

  cat(m, "rows:", nrow(df_target[df_target$mut == m, ]), "\n")

}


library(dplyr)

library(tidyr)

library(stringr)

library(glmmTMB)

library(emmeans)

library(ggplot2)

library(patchwork)

library(showtext)

library(sysfonts)


project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)

input_dir <- file.path(project_root, "inputs")

out_dir <- project_root


mutation_file <- file.path(input_dir, "all_pass_snps_per_sample.tsv")

callable_file <- file.path(input_dir, "callable_base_counts.tsv")


raw_mut_levels <- c(

  "A>C", "A>G", "A>T",

  "C>A", "C>G", "C>T",

  "G>A", "G>C", "G>T",

  "T>A", "T>C", "T>G"

)


mut_levels <- c(

  "C>T/G>A",

  "C>A/G>T",

  "C>G/G>C",

  "T>A/A>T",

  "T>C/A>G",

  "T>G/A>C"

)


mut_label <- paste(mut_levels, collapse = " + ")


mut_plot_labels <- c(

  "C>T/G>A" = "C>T",

  "C>A/G>T" = "C>A",

  "C>G/G>C" = "C>G",

  "T>A/A>T" = "T>A",

  "T>C/A>G" = "T>C",

  "T>G/A>C" = "T>G"

)


dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


cond_levels <- c("H", "L")


base_family <- "Times New Roman"


font_add(

  family = "Times New Roman",

  regular = "/usr/share/fonts/TTF/Times.TTF",

  bold    = "/usr/share/fonts/TTF/Timesbd.TTF",

  italic  = "/usr/share/fonts/TTF/Timesi.TTF",

  bolditalic = "/usr/share/fonts/TTF/Timesbi.TTF"

)


showtext_auto()


pval_to_stars <- function(p) {

  case_when(

    is.na(p) ~ "",

    p < 0.001 ~ "***",

    p < 0.01  ~ "**",

    p < 0.05  ~ "*",

    TRUE ~ ""

  )

}


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


df <- read.delim(

  mutation_file,

  sep = "\t",

  header = TRUE,

  stringsAsFactors = FALSE

)


callable <- read.delim(

  callable_file,

  sep = "\t",

  header = TRUE,

  stringsAsFactors = FALSE

)


required_mut_cols <- c("sample", "REF", "ALT")

missing_mut_cols <- setdiff(required_mut_cols, names(df))


if (length(missing_mut_cols) > 0) {

  stop(

    "Missing required columns in mutation file: ",

    paste(missing_mut_cols, collapse = ", ")

  )

}


required_callable_cols <- c("sample", "A", "C", "G", "T")

missing_callable_cols <- setdiff(required_callable_cols, names(callable))


if (length(missing_callable_cols) > 0) {

  stop(

    "Missing required columns in callable file: ",

    paste(missing_callable_cols, collapse = ", ")

  )

}


if (anyDuplicated(callable$sample) > 0) {

  stop("Duplicate sample names found in callable file.")

}


mutation_samples <- unique(df$sample)

callable_samples <- unique(callable$sample)


samples_to_use <- intersect(mutation_samples, callable_samples)


if (length(samples_to_use) == 0) {

  stop("No overlapping samples between mutation file and callable file.")

}


extra_callable_samples <- setdiff(callable_samples, mutation_samples)


if (length(extra_callable_samples) > 0) {

  message(

    "Ignoring callable-only samples not present in mutation file: ",

    paste(extra_callable_samples, collapse = ", ")

  )

}


missing_callable_samples <- setdiff(mutation_samples, callable_samples)


if (length(missing_callable_samples) > 0) {

  stop(

    "These mutation samples are missing from callable file: ",

    paste(missing_callable_samples, collapse = ", ")

  )

}


df <- df %>%

  filter(sample %in% samples_to_use)


callable <- callable %>%

  filter(sample %in% samples_to_use)


sample_meta <- data.frame(

  sample = samples_to_use,

  stringsAsFactors = FALSE

) %>%

  mutate(

    rep_id = str_match(sample, "^.+_([0-9]+)_([HL])$")[, 2],

    cond   = str_match(sample, "^.+_([0-9]+)_([HL])$")[, 3]

  ) %>%

  filter(

    !is.na(rep_id),

    !is.na(cond),

    cond %in% cond_levels

  ) %>%

  mutate(

    rep_id = factor(rep_id),

    cond = factor(cond, levels = cond_levels)

  )


if (nrow(sample_meta) == 0) {

  stop("No samples matched expected pattern like WRS_1_H / WRS_1_L.")

}


df2 <- df %>%

  mutate(

    REF = toupper(REF),

    ALT = toupper(ALT),


    is_snv = REF %in% c("A", "C", "G", "T") &

      ALT %in% c("A", "C", "G", "T"),


    raw_mut = paste0(REF, ">", ALT),

    mut_type = collapse_mut(raw_mut)

  ) %>%

  filter(

    is_snv,

    raw_mut %in% raw_mut_levels,

    mut_type %in% mut_levels

  ) %>%

  left_join(sample_meta, by = "sample") %>%

  filter(

    !is.na(rep_id),

    !is.na(cond)

  ) %>%

  mutate(

    mut_type = factor(mut_type, levels = mut_levels)

  )


mut_counts <- df2 %>%

  count(sample, rep_id, cond, mut_type, name = "mut_n")


mut_counts <- sample_meta %>%

  tidyr::expand_grid(

    mut_type = factor(mut_levels, levels = mut_levels)

  ) %>%

  left_join(

    mut_counts,

    by = c("sample", "rep_id", "cond", "mut_type")

  ) %>%

  mutate(

    mut_n = ifelse(is.na(mut_n), 0L, mut_n),

    mut_n = as.integer(mut_n)

  )


df_model <- mut_counts %>%

  left_join(callable, by = "sample") %>%

  mutate(

    A = as.numeric(A),

    C = as.numeric(C),

    G = as.numeric(G),

    T = as.numeric(T),


    callable_n = case_when(

      as.character(mut_type) == "C>A/G>T" ~ C + G,

      as.character(mut_type) == "C>G/G>C" ~ C + G,

      as.character(mut_type) == "C>T/G>A" ~ C + G,

      as.character(mut_type) == "T>A/A>T" ~ T + A,

      as.character(mut_type) == "T>C/A>G" ~ T + A,

      as.character(mut_type) == "T>G/A>C" ~ T + A,

      TRUE ~ NA_real_

    ),


    cond = factor(cond, levels = cond_levels),

    mut_type = factor(mut_type, levels = mut_levels),

    rep_id = factor(rep_id),


    observed_VR = mut_n / callable_n

  ) %>%

  filter(

    !is.na(callable_n),

    callable_n > 0,

    !is.na(cond),

    !is.na(mut_type),

    !is.na(rep_id)

  )


if (nrow(df_model) == 0) {

  stop("No rows left for model after target-substitution filtering: ", mut_label)

}


fit <- glmmTMB(

  mut_n ~ cond * mut_type + offset(log(callable_n)) + (1 | rep_id),

  data = df_model,

  family = nbinom2,

  control = glmmTMBControl(

    optCtrl = list(iter.max = 1e5, eval.max = 1e5)

  )

)


emm_vr <- emmeans(

  fit,

  ~ cond | mut_type,

  type = "response",

  offset = 0

)


emm_vr_dt <- as.data.frame(summary(emm_vr))


if (!"response" %in% names(emm_vr_dt)) {

  if ("rate" %in% names(emm_vr_dt)) {

    emm_vr_dt <- emm_vr_dt %>% rename(response = rate)

  } else {

    stop("Could not find response/rate column in emmeans output. Check names(emm_vr_dt).")

  }

}


if (!all(c("asymp.LCL", "asymp.UCL") %in% names(emm_vr_dt))) {

  if (all(c("lower.CL", "upper.CL") %in% names(emm_vr_dt))) {

    emm_vr_dt <- emm_vr_dt %>%

      rename(

        asymp.LCL = lower.CL,

        asymp.UCL = upper.CL

      )

  } else {

    stop("Could not find confidence interval columns in emmeans output.")

  }

}


emm_vr_dt <- emm_vr_dt %>%

  mutate(

    mut_type = factor(mut_type, levels = mut_levels),

    cond = factor(cond, levels = cond_levels)

  ) %>%

  arrange(mut_type, cond)


contrast_vr <- contrast(

  emm_vr,

  method = "revpairwise",

  by = "mut_type",

  type = "response"

)


posthoc_L_vs_H <- as.data.frame(summary(

  contrast_vr,

  infer = TRUE,

  adjust = "none"

))


if (!"ratio" %in% names(posthoc_L_vs_H)) {

  if ("rate.ratio" %in% names(posthoc_L_vs_H)) {

    posthoc_L_vs_H <- posthoc_L_vs_H %>%

      rename(ratio = rate.ratio)

  } else if ("response" %in% names(posthoc_L_vs_H)) {

    posthoc_L_vs_H <- posthoc_L_vs_H %>%

      rename(ratio = response)

  } else {

    stop("Could not find fold-change column. Check names(posthoc_L_vs_H).")

  }

}


if (!all(c("asymp.LCL", "asymp.UCL") %in% names(posthoc_L_vs_H))) {

  if (all(c("lower.CL", "upper.CL") %in% names(posthoc_L_vs_H))) {

    posthoc_L_vs_H <- posthoc_L_vs_H %>%

      rename(

        asymp.LCL = lower.CL,

        asymp.UCL = upper.CL

      )

  } else {

    stop("Could not find confidence interval columns in posthoc output.")

  }

}


posthoc_L_vs_H <- posthoc_L_vs_H %>%

  mutate(

    mut_type = factor(mut_type, levels = mut_levels)

  ) %>%

  arrange(mut_type) %>%

  rename(p_value = p.value) %>%

  mutate(

    p_holm = p.adjust(p_value, method = "holm"),

    p_label = pval_to_stars(p_holm)

  ) %>%

  select(

    any_of(c(

      "contrast", "mut_type", "ratio", "SE", "df",

      "asymp.LCL", "asymp.UCL", "null", "z.ratio", "t.ratio",

      "p_value", "p_holm", "p_label"

    )),

    everything()

  )


sig_dt <- posthoc_L_vs_H %>%

  select(mut_type, p_value, p_holm, p_label)


y_range_A <- range(

  c(emm_vr_dt$asymp.LCL, emm_vr_dt$asymp.UCL),

  na.rm = TRUE

)


y_offset_A <- diff(y_range_A) * 0.04


sig_pos_A <- emm_vr_dt %>%

  group_by(mut_type) %>%

  summarise(

    y_pos = max(asymp.UCL, na.rm = TRUE) + y_offset_A,

    .groups = "drop"

  ) %>%

  left_join(sig_dt, by = "mut_type") %>%

  mutate(

    mut_type = factor(mut_type, levels = mut_levels)

  )


p_est_vr_bar <- ggplot(

  emm_vr_dt,

  aes(x = mut_type, y = response, fill = cond)

) +

  geom_col(

    position = position_dodge(width = 0.8),

    width = 0.7,

    color = "black",

    linewidth = 0.2

  ) +

  geom_errorbar(

    aes(ymin = asymp.LCL, ymax = asymp.UCL),

    position = position_dodge(width = 0.8),

    width = 0.2,

    linewidth = 0.4

  ) +

  geom_text(

    data = sig_pos_A,

    aes(x = mut_type, y = y_pos, label = p_label),

    inherit.aes = FALSE,

    size = 5,

    family = base_family

  ) +

  scale_x_discrete(labels = mut_plot_labels, drop = FALSE) +

  scale_y_continuous(

    expand = expansion(mult = c(0, 0.10))

  ) +

  theme_classic(base_size = 15, base_family = base_family) +

  labs(

    x = "Substitution",

    y = "Estimated VR",

    fill = "Condition"

  ) +

  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))


posthoc_fc_dt <- posthoc_L_vs_H %>%

  mutate(

    mut_type = factor(mut_type, levels = mut_levels)

  )


p_fc_vr <- ggplot(

  posthoc_fc_dt,

  aes(x = mut_type, y = ratio)

) +

  geom_hline(yintercept = 1, linetype = "dashed") +

  geom_col(

    width = 0.7,

    color = "black",

    linewidth = 0.2

  ) +

  geom_errorbar(

    aes(ymin = asymp.LCL, ymax = asymp.UCL),

    width = 0.2,

    linewidth = 0.4

  ) +

  geom_text(

    aes(

      y = asymp.UCL * 1.05,

      label = p_label

    ),

    size = 5,

    family = base_family

  ) +

  scale_x_discrete(labels = mut_plot_labels, drop = FALSE) +

  scale_y_continuous(

    limits = c(0, NA),

    expand = expansion(mult = c(0, 0.18))

  ) +

  theme_classic(base_size = 15, base_family = base_family) +

  labs(

    x = "Substitution",

    y = "Estimated VR fold change\nLow DO / High DO"

  ) +

  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))


p_est_vr_bar_clean <- p_est_vr_bar +

  labs(tag = "A") +

  theme(

    legend.position = "right",

    plot.tag = element_text(size = 18, face = "bold", family = base_family),

    plot.tag.position = c(-0.03, 1.03),

    plot.margin = margin(t = 12, r = 8, b = 5, l = 18)

  )


p_fc_vr_clean <- p_fc_vr +

  labs(tag = "B") +

  theme(

    legend.position = "none",

    plot.tag = element_text(size = 18, face = "bold", family = base_family),

    plot.tag.position = c(-0.03, 1.03),

    plot.margin = margin(t = 12, r = 8, b = 5, l = 18)

  )


p_combined_vr <- p_est_vr_bar_clean + p_fc_vr_clean +

  plot_layout(widths = c(1.25, 1))


ggsave(

  file.path(out_dir, "nbinom_estimated_VR_and_fold_change.pdf"),

  p_combined_vr,

  width = 13,

  height = 5

)


cat("Done.\n")

cat("Results saved to:", out_dir, "\n")

cat("Samples used:", paste(samples_to_use, collapse = ", "), "\n")

cat("Target substitutions:", mut_label, "\n")

cat("Rows used in model:", nrow(df_model), "\n")

for (m in mut_levels) {

  cat(m, "rows:", nrow(df_model[df_model$mut_type == m, ]), "\n")

}

for (m in mut_levels) {

  cat("Total", m, "mutations:", sum(df_model$mut_n[df_model$mut_type == m]), "\n")

}
