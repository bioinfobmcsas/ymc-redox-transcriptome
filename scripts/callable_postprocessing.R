library(data.table)
library(stringr)

project_root <- normalizePath(if (dir.exists("inputs")) "." else "..", mustWork = TRUE)
input_dir <- file.path(project_root, "inputs")
dir_path <- file.path(input_dir, "loci_counts")

cond_ref <- "H"
cond_test <- "L"

files <- list.files(
  dir_path,
  pattern = "\\.csv\\.gz$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No .csv.gz files found in: ", dir_path)
}

read_one <- function(f) {
  dt <- fread(f)

  dt[, sample := sub("\\.ref_alt_counts\\.csv\\.gz$", "", basename(f))]
  dt[, sample := sub("\\..*$", "", sample)]
  dt[, coverage := COUNT_REF + COUNT_ALT]

  dt
}

dt_all <- rbindlist(lapply(files, read_one), use.names = TRUE, fill = TRUE)

required_cols <- c("CHR", "POS", "REF", "ALT", "COUNT_REF", "COUNT_ALT")
missing_cols <- setdiff(required_cols, names(dt_all))

if (length(missing_cols) > 0) {
  stop("Missing columns in loci-count files: ", paste(missing_cols, collapse = ", "))
}

dt_all[, mut := paste0(REF, ">", ALT)]
dt_all[, cond := str_sub(sample, -1)]
dt_all[, cond := factor(cond, levels = c(cond_ref, cond_test))]
dt_all[, id := paste0(CHR, ":", POS, ":", REF, ":", ALT)]

comp <- c(A = "T", T = "A", C = "G", G = "C")

dt_all[, mut_type := fifelse(
  REF %in% c("C", "T"),
  paste0(REF, ">", ALT),
  paste0(comp[REF], ">", comp[ALT])
)]

mut_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
dt_all[, mut_type := factor(mut_type, levels = mut_levels)]

dt_filt <- copy(dt_all)

dt_sum <- dt_filt[
  ,
  .(
    COUNT_REF = sum(COUNT_REF, na.rm = TRUE),
    COUNT_ALT = sum(COUNT_ALT, na.rm = TRUE)
  ),
  by = .(mut_type, id, cond)
]

dt_sum[, total := COUNT_REF + COUNT_ALT]
dt_sum <- dt_sum[total > 0]
dt_sum[, frac_alt := COUNT_ALT / total]

dt_wide <- dcast(
  dt_sum,
  mut_type + id ~ cond,
  value.var = "frac_alt"
)

dt_wide[, delta := get(cond_test) - get(cond_ref)]
dt_delta <- dt_wide[!is.na(delta)]

test_one_group <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) < 2) {
    return(data.table(
      n_sites = length(x),
      mean_delta = mean(x),
      median_delta = median(x),
      sd_delta = sd(x),
      t_statistic = NA_real_,
      t_p_value = NA_real_,
      t_ci_low = NA_real_,
      t_ci_high = NA_real_
    ))
  }

  tt <- t.test(x, mu = 0)

  data.table(
    n_sites = length(x),
    mean_delta = mean(x),
    median_delta = median(x),
    sd_delta = sd(x),
    t_statistic = unname(tt$statistic),
    t_p_value = tt$p.value,
    t_ci_low = tt$conf.int[1],
    t_ci_high = tt$conf.int[2]
  )
}

res <- dt_delta[
  ,
  test_one_group(delta),
  by = mut_type
]

dt_filt[, replicate := str_match(sample, "^WRS_([0-9]+)_[HL]$")[, 2]]

dt_sample_rate <- dt_filt[
  ,
  .(
    COUNT_REF = sum(COUNT_REF, na.rm = TRUE),
    COUNT_ALT = sum(COUNT_ALT, na.rm = TRUE)
  ),
  by = .(sample, replicate, cond, mut_type)
]

dt_sample_rate[, total := COUNT_REF + COUNT_ALT]
dt_sample_rate <- dt_sample_rate[total > 0]
dt_sample_rate[, alt_fraction := COUNT_ALT / total]

fwrite(dt_sum, file.path(input_dir, "loci_counts_summed_by_site_cond.csv.gz"))
fwrite(dt_delta, file.path(input_dir, "loci_delta_by_site.csv.gz"))
fwrite(res, file.path(input_dir, "loci_delta.csv"))
fwrite(dt_sample_rate, file.path(input_dir, "loci_per_sample_alt_fraction.csv"))
