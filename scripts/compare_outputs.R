# Compare outputs of the original workbook and the cleaned pipeline.
#
# 1. Re-run the ORIGINAL workbook with the current input files, so both versions
#    see identical inputs (the old outputs on disk may predate later edits to the
#    flag files).
# 2. Run fluorescence_correction_pipeline.Rmd.
# 3. Set the file names below and run this script.
#
# Expected differences (everything else should match exactly):
#   - New column `profile_correction` (only in new).
#   - Saturated profiles (2022-03-09, 2022-09-09 KC10): `ctdNum` is now filled in
#     (was NA); duplicate rows from replicate casts are gone; profiles at other
#     stations on those dates are no longer overwritten with NA.
#   - Size-fractionated sums: the old code de-duplicated on the `sum` value
#     across all samples, which would drop a sample whose sum happened to equal
#     another's. Probably no effect, but if you see differences in `sum` or
#     `chl_comb`, this is why.
#   - If quench_correction_fits.csv (old) and .xlsx (new) differ, NPQ profiles
#     will differ.

library(dplyr)
library(readr)
library(here)

old_final <- here("outputs", "flu_qc_2025-11-25.csv")
new_final <- here("outputs", "flu_qc_2026-XX-XX.csv")  # set to the new file
old_cal   <- here("outputs", "ctd_fluorescence_calibrations_stats_2025-11-18.csv")
new_cal   <- here("outputs", "ctd_fluorescence_calibration_stats_2026-XX-XX.csv")

tol <- 1e-8

read_output <- function(path) {
  read_csv(path, guess_max = 1e6, show_col_types = FALSE) %>%
    select(-any_of(c("...1", "X")))  # row-name column from write.csv()
}

differs <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    out <- abs(a - b) > tol
  } else {
    a <- as.character(a)
    b <- as.character(b)
    out <- a != b
  }
  na_out <- is.na(out)
  out[na_out] <- is.na(a[na_out]) != is.na(b[na_out])
  out
}

compare_tables <- function(old, new, keys, label) {
  cat("\n==========", label, "==========\n")
  cat("Rows: old", nrow(old), "| new", nrow(new), "\n")

  only_old_cols <- setdiff(names(old), names(new))
  only_new_cols <- setdiff(names(new), names(old))
  if (length(only_old_cols)) cat("Columns only in old:", only_old_cols, "\n")
  if (length(only_new_cols)) cat("Columns only in new:", only_new_cols, "\n")

  # Number duplicate keys so duplicated rows line up one-to-one
  add_dup <- function(d) {
    d %>% group_by(across(all_of(keys))) %>% mutate(.dup = row_number()) %>% ungroup()
  }
  joined <- full_join(add_dup(old) %>% mutate(.in_old = TRUE),
                      add_dup(new) %>% mutate(.in_new = TRUE),
                      by = c(keys, ".dup"), suffix = c(".old", ".new"))

  only_old <- filter(joined, is.na(.in_new))
  only_new <- filter(joined, is.na(.in_old))
  if (nrow(only_old)) {
    cat("\nRows only in old:", nrow(only_old), "\n")
    print(count(only_old, across(all_of(keys[1]))), n = 50)
  }
  if (nrow(only_new)) {
    cat("\nRows only in new:", nrow(only_new), "\n")
    print(count(only_new, across(all_of(keys[1]))), n = 50)
  }

  both <- filter(joined, .in_old, .in_new)
  common <- setdiff(intersect(names(old), names(new)), keys)
  mismatch <- lapply(common, function(col) {
    d <- differs(both[[paste0(col, ".old")]], both[[paste0(col, ".new")]])
    tibble(column = col, n_diff = sum(d), key = both[[keys[1]]][d])
  })
  summary_tbl <- tibble(column = common,
                        n_diff = vapply(mismatch, function(m) if (nrow(m)) m$n_diff[1] else 0L, integer(1)),
                        example_keys = vapply(mismatch, function(m) paste(head(unique(m$key), 5), collapse = ", "), ""))

  cat("\nMatched rows:", nrow(both), "\n")
  if (all(summary_tbl$n_diff == 0)) {
    cat("All shared columns identical (tolerance", tol, ").\n")
  } else {
    cat("Columns with differences:\n")
    print(filter(summary_tbl, n_diff > 0), n = 50, width = Inf)
  }
  invisible(list(joined = joined, summary = summary_tbl))
}

res_cal   <- compare_tables(read_output(old_cal), read_output(new_cal),
                            keys = c("ctdNum", "calibration_type"), label = "Calibration stats")
res_final <- compare_tables(read_output(old_final), read_output(new_final),
                            keys = c("id", "pres"), label = "Final fluorescence product")

# To inspect a specific profile:
# res_final$joined %>% filter(id == "2022-03-09-KC10-18032") %>% View()
