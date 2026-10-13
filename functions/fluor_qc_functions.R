# Helper functions for the CTD fluorescence correction pipeline
# (fluorescence_correction_pipeline.Rmd). Sourced at the top of the notebook.

# ---- Manual decisions --------------------------------------------------------

# Allowed decision types in manual_decisions.csv, and what `id` means for each.
decision_types <- c(
  exclude_fit_point           = "Cast ID (date-station-ctdNum) + pres: drop one matchup from the cast-specific fits",
  exclude_single_depth_cast   = "Cast ID: drop a 1-2 sample cast from the sensor-wide calibration",
  exclude_single_depth_sensor = "Sensor serial: drop all 1-2 sample casts from that sensor",
  separate_calibration        = "Cast ID: fit these casts as their own calibration group for that sensor",
  replace_saturated           = "date-station: replace the profile with interpolated discrete chlorophyll"
)

read_decisions <- function(path) {
  d <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()))
  unknown <- setdiff(unique(d$decision), names(decision_types))
  if (length(unknown) > 0) {
    stop("Unknown decision type(s) in ", path, ": ", paste(unknown, collapse = ", "))
  }
  d
}

decision_ids <- function(decisions, type) {
  decisions$id[decisions$decision == type]
}

# Warn when a manual decision doesn't match anything in the data
# (catches typos and decisions that went stale after a data re-download).
warn_unmatched <- function(ids, available, label) {
  missing <- setdiff(ids, available)
  if (length(missing) > 0) {
    warning(label, ": not found in data: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(missing)
}

# Which calibration group each row belongs to. Most sensors have one ("standard").
# Sensors with casts listed under `separate_calibration` get two groups: those
# casts ("april_only") and everything else ("excl_april"). The names are kept
# from the original notebook so outputs stay comparable.
assign_calibration_type <- function(df, sep_ids) {
  sep_sensors <- unique(sub(".*-", "", sep_ids))
  df %>%
    dplyr::mutate(calibration_type = dplyr::case_when(
      !ctdNum %in% sep_sensors ~ "standard",
      id %in% sep_ids          ~ "april_only",
      TRUE                     ~ "excl_april"
    ))
}

# ---- Regression helpers ------------------------------------------------------

# Overall linear fit chl_comb ~ f_dm (unrounded statistics).
fit_line <- function(d) {
  m <- lm(chl_comb ~ f_dm, data = d)
  s <- summary(m)
  tibble::tibble(
    slope     = unname(coef(m)[2]),
    intercept = unname(coef(m)[1]),
    r2        = s$r.squared,
    adj_r2    = s$adj.r.squared,
    p_value   = s$coefficients[2, 4],
    n         = length(m$residuals)
  )
}

# Rounded fit used for cast-by-cast models; NA if fewer than 3 points.
fit_cast_once <- function(d) {
  if (nrow(d) < 3) {
    return(list(intercept = NA, slope = NA, r2 = NA, p.value = NA))
  }
  m <- lm(chl_comb ~ f_dm, data = d)
  s <- summary(m)
  list(
    intercept = round(unname(coef(m)[1]), 2),
    slope     = round(unname(coef(m)[2]), 2),
    r2        = round(s$adj.r.squared, 2),
    p.value   = round(s$coefficients[2, 4], 5)
  )
}

# Fit one cast with and without the 1 m sample, and keep the fit with the lower
# (rounded) p-value. Dropping the surface handles casts where NPQ quenching makes
# the surface CTD value artificially low. Note this compares p-values from fits
# with different n, so it is a heuristic, not a formal model comparison.
select_cast_model <- function(d) {
  no_surf   <- fit_cast_once(dplyr::filter(d, pres != 1))
  with_surf <- fit_cast_once(d)
  n_no   <- sum(d$pres != 1)
  n_with <- nrow(d)
  can_no   <- n_no >= 3
  can_with <- n_with >= 3

  if (!can_no && !can_with) {
    use_surface <- NA
    best <- list(intercept = NA, slope = NA, r2 = NA, p.value = NA)
    nobs <- NA
  } else if (!can_no) {
    use_surface <- TRUE
    best <- with_surf
    nobs <- n_with
  } else if (!can_with) {
    use_surface <- FALSE
    best <- no_surf
    nobs <- n_no
  } else {
    use_surface <- with_surf$p.value < no_surf$p.value
    best <- purrr::map2(with_surf, no_surf, ~ ifelse(use_surface, .x, .y))
    nobs <- ifelse(use_surface, n_with, n_no)
  }

  surface_available <- any(d$pres == 1)

  tibble::tibble(
    intercept   = best$intercept,
    slope       = best$slope,
    r2          = best$r2,
    p.value     = best$p.value,
    nobs        = nobs,
    use_surface = use_surface,
    quenched    = surface_available && !isTRUE(use_surface)
  )
}

# ---- NPQ hybrid correction ---------------------------------------------------

# For one profile: fit a loess through the discrete samples (chl_comb ~ pres),
# then use the loess above the "transition depth" and the slope-corrected
# fluorescence below it. The transition depth is the shallowest depth below 5 m
# where the two agree to within `threshold` (mg m-3); if they never agree, the
# loess is used for the whole profile. Below the deepest discrete sample the
# loess is NA, and the pipeline falls back to f_calibrated there.
npq_hybrid <- function(profile, threshold, span) {
  discrete <- dplyr::filter(profile, !is.na(chl_comb))
  if (nrow(discrete) < 3) {
    return(dplyr::mutate(profile,
                         chl_fitted_raw = NA_real_,
                         transition_depth = NA_real_,
                         f_npq = f_calibrated))
  }
  fit <- loess(chl_comb ~ pres, data = discrete, span = span)
  profile$chl_fitted_raw <- predict(fit, newdata = profile)
  profile %>%
    dplyr::mutate(
      npq_diff = abs(chl_fitted_raw - f_calibrated),
      transition_depth = {
        idx <- which(npq_diff < threshold & pres > 5)
        if (length(idx) > 0) min(pres[idx]) else max(pres)
      },
      f_npq = dplyr::if_else(pres < transition_depth, chl_fitted_raw, f_calibrated)
    )
}

# ---- Saturated-profile interpolation -----------------------------------------

# Loess-interpolate discrete chlorophyll profiles to 1 m resolution.
# `discrete` needs columns date, station, pres, chla.
interp_discrete_profiles <- function(discrete, span) {
  wide <- discrete %>%
    dplyr::mutate(key = paste(date, station, sep = "_")) %>%
    dplyr::select(key, pres, chla) %>%
    tidyr::pivot_wider(names_from = key, values_from = chla)

  depth <- wide$pres
  x_new <- seq(1, max(depth), by = 1)
  interp <- rioja::interp.dataset(y = dplyr::select(wide, -pres), x = depth,
                                  xout = x_new, method = "loess", span = span)

  as.data.frame(interp) %>%
    dplyr::mutate(pres = x_new) %>%
    tidyr::pivot_longer(-pres, names_to = "key", values_to = "f_interpolated") %>%
    tidyr::separate(key, into = c("date", "station"), sep = "_") %>%
    dplyr::mutate(date = as.Date(date),
                  f_interpolated = round(f_interpolated, 2))
}

# ---- Diagnostic plots --------------------------------------------------------

integrate_trap <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

# Profile + scatter for one cast, to diagnose poor fluorescence/discrete matchups
# (often depth mismatches around sharp peaks or an SCM).
# cast_id: "date-station-ctdNum" (or "date-station" for all sensors that day).
plot_cast_detailed <- function(cast_id, f_dm, f_fit_multi, cast_models, sensor_medians) {
  parts <- strsplit(cast_id, "-")[[1]]
  cast_date <- as.Date(paste(parts[1:3], collapse = "-"))
  cast_station <- parts[4]

  pick_cast <- function(df) {
    df <- dplyr::filter(df, date == cast_date, station == cast_station)
    if (length(parts) == 5) df <- dplyr::filter(df, ctdNum == parts[5])
    df
  }
  profile_data <- pick_cast(f_dm)
  cast_data <- pick_cast(f_fit_multi)
  if (nrow(profile_data) == 0) stop("No profile data found for cast: ", cast_id)

  sensor <- profile_data$ctdNum[1]
  sensor_med <- dplyr::filter(sensor_medians, ctdNum == sensor)

  cast_info <- dplyr::filter(cast_models, id == cast_id)
  if (nrow(cast_info) == 0) {
    warning("Cast ", cast_id, " not in cast_models; using surface by default.")
    use_surface <- TRUE
  } else {
    use_surface <- isTRUE(cast_info$use_surface[1])
  }

  fit_data <- cast_data %>% dplyr::filter(!is.na(chl_comb), !is.na(f_dm))
  if (!use_surface) fit_data <- dplyr::filter(fit_data, pres != 1)

  m <- lm(chl_comb ~ f_dm, data = fit_data)
  eq_label <- sprintf(
    "Cast: y = %.2fx + %.2f\nMedian: y = %.2fx + %.2f\nR² = %.2f, p = %.4f",
    coef(m)[2], coef(m)[1], sensor_med$median_slope, sensor_med$median_intercept,
    summary(m)$r.squared, summary(m)$coefficients[2, 4]
  )

  profile_plot <- profile_data %>%
    dplyr::filter(pres < 50) %>%
    ggplot2::ggplot(ggplot2::aes(y = pres)) +
    ggplot2::geom_ribbon(ggplot2::aes(xmin = f_lower, xmax = f_upper),
                         fill = "#2E86AB", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(x = f_dm, color = "CTD Fluorescence"),
                       linewidth = 1, orientation = "y") +
    ggplot2::geom_point(data = dplyr::filter(cast_data, pres < 50),
                        ggplot2::aes(x = chl_comb, color = "Discrete Samples"), size = 3) +
    ggplot2::scale_y_reverse(name = "Pressure (dbar)") +
    ggplot2::scale_x_continuous(name = "Chlorophyll (µg/L)") +
    ggplot2::scale_color_manual(name = "", values = c("CTD Fluorescence" = "#2E86AB",
                                                      "Discrete Samples" = "#A23B72")) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = "A) Fluorescence Profile")

  lims <- range(c(cast_data$f_dm, cast_data$chl_comb), na.rm = TRUE)
  buffer <- diff(lims) * 0.05
  axis_limits <- c(max(0, lims[1] - buffer), lims[2] + buffer)

  scatter_plot <- cast_data %>%
    dplyr::filter(!is.na(chl_comb), !is.na(f_dm)) %>%
    ggplot2::ggplot(ggplot2::aes(x = f_dm, y = chl_comb, color = as.factor(pres))) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    ggplot2::geom_abline(slope = sensor_med$median_slope,
                         intercept = sensor_med$median_intercept, color = "blue", linewidth = 1.2) +
    ggplot2::geom_smooth(data = fit_data, ggplot2::aes(x = f_dm, y = chl_comb),
                         method = "lm", color = "red", fill = "red", alpha = 0.2,
                         inherit.aes = FALSE) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_brewer(name = "Pressure\n(dbar)", palette = "Set2") +
    ggplot2::scale_x_continuous(name = "CTD Fluorescence (µg/L)", limits = axis_limits) +
    ggplot2::scale_y_continuous(name = "Discrete Chlorophyll (µg/L)", limits = axis_limits) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = "B) CTD vs Discrete Chlorophyll",
                  subtitle = paste0("Red = cast fit | Blue = sensor median | Dashed = 1:1",
                                    if (!use_surface) " | Surface excluded" else "")) +
    ggplot2::annotate("text", x = axis_limits[1], y = axis_limits[2], label = eq_label,
                      hjust = -0.05, vjust = 1.1, size = 3.5, fontface = "italic")

  profile_plot + scatter_plot +
    patchwork::plot_annotation(title = paste0("Cast: ", cast_id, " | Sensor: ", sensor))
}

# Show the NPQ hybrid correction for one profile, with 0-30 m integrals.
# Use this to choose `threshold` and `span` for quench_correction_fits.xlsx.
plot_npq_correction <- function(profile_id, data, threshold = 0.1, span = 0.75,
                                max_depth = 50) {
  prof <- npq_hybrid(dplyr::filter(data, id == profile_id), threshold, span)
  trans <- unique(prof$transition_depth)

  integ <- prof %>% dplyr::filter(pres >= 0, pres <= 30) %>% dplyr::arrange(pres)
  x_lab <- max(prof$f_calibrated, na.rm = TRUE) * 0.95
  ann <- tibble::tibble(
    x = x_lab, y = c(28, 32),
    label = c(paste0("Slope: ", round(integrate_trap(integ$pres, integ$f_calibrated), 1), " mg/m²"),
              paste0("Hybrid: ", round(integrate_trap(integ$pres, integ$f_npq), 1), " mg/m²")),
    color = c("Slope-Corrected", "Hybrid Correction")
  )

  prof %>%
    dplyr::filter(pres < max_depth) %>%
    ggplot2::ggplot(ggplot2::aes(y = pres)) +
    ggplot2::geom_line(ggplot2::aes(x = f_dm, color = "Raw CTD Fluorescence"), orientation = "y") +
    ggplot2::geom_line(ggplot2::aes(x = f_calibrated, color = "Slope-Corrected"), orientation = "y") +
    ggplot2::geom_line(ggplot2::aes(x = f_npq, color = "Hybrid Correction"),
                       linewidth = 1.2, orientation = "y") +
    ggplot2::geom_point(data = dplyr::filter(prof, pres < max_depth, !is.na(chl_comb)),
                        ggplot2::aes(x = chl_comb, color = "Discrete Samples"), size = 3) +
    ggplot2::geom_text(data = ann, ggplot2::aes(x = x, y = y, label = label, color = color),
                       hjust = 1, size = 3.5, fontface = "bold", show.legend = FALSE) +
    ggplot2::scale_y_reverse(name = "Pressure (dbar)") +
    ggplot2::scale_x_continuous(name = "Chlorophyll (µg/L)") +
    ggplot2::scale_color_manual(name = "", values = c(
      "Raw CTD Fluorescence" = "#2E86AB", "Slope-Corrected" = "#F18F01",
      "Hybrid Correction" = "#6A4C93", "Discrete Samples" = "#A23B72")) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = paste("NPQ correction:", profile_id),
                  subtitle = if (!all(is.na(trans))) {
                    paste("Transition depth:", round(trans[1], 1), "dbar | span:", span,
                          "| threshold:", threshold)
                  })
}
