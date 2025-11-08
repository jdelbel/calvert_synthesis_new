# ==============================================================================
# CTD Fluorescence Cast-by-Cast QC Tool
# ==============================================================================
# This script loads CTD fluorescence data, discrete chlorophyll samples,
# calculates cast-by-cast calibrations, and provides an interactive QC tool
# to review and flag problematic casts.
#
# Usage: Source this script and run browse_casts()
# ==============================================================================

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)
library(readr)
library(lubridate)

# ==============================================================================
# STEP 1: Load and prepare data
# ==============================================================================

cat("Loading data...\n")

# Downloading baseline corrected chlorophyll fluorescence profiles
f <- read_csv(here("files", "8_binAvg-1762199720021.csv"))

# Flagging file
flags <- read_csv(here("outputs","fluorescence_QC_flags_20251103.csv"))

# Downloading chlorophyll data for joining
c <- read_csv(here("files", "2025-11-03_HakaiData_chlorophyll.csv"))

cat("Data loaded successfully.\n")

# ==============================================================================
# STEP 2: Wrangle CTD profiles
# ==============================================================================

cat("Wrangling CTD profiles...\n")

# Wrangling CTD profiles, setting date column and renaming columns 
f <- f %>%
  filter(`Cast Direction Flag` == "d") %>%  #downcast data
  mutate(date = lubridate::date(`Measurement time`)) %>%
  mutate(year = lubridate::year(`Measurement time`)) %>%
  select(castpk = `Cast PK`,
         hakai_id = `Hakai ID`,
         Cruise,
         ctdNum = `CTD serial number`,
         station = Station,
         lat = Latitude...11,
         long = Longitude...12,
         time = `Measurement time`,
         date,
         year,
         pres = `Pressure (dbar)`,
         flu = `Fluorometry Chlorophyll (ug/L)`,
         flu_flag = `Fluorometry Chlorophyll flag`) %>% 
  mutate(time = lubridate::with_tz(time, tzone = "America/Vancouver")) 

# Removing bad or suspect profiles, but leaving shallow casts.
f <- f %>% 
  left_join(flags) %>% 
  filter(flag == "AV") %>% 
  filter(station %in% c("KC10", "FZH01"))

# For now, doing a daily mean on fluorescence profiles
f_dm <- f %>% 
  group_by(date, station, pres, ctdNum) %>% 
  summarise(f_dm = mean(flu),
            n_prof = n()) %>% 
  ungroup() %>% 
  mutate(f_dm = round(f_dm, 2)) %>% 
  unite(id, c(date, station), sep = "-", remove = F)

cat("CTD profiles prepared.\n")


# ==============================================================================
# STEP 3: Process discrete chlorophyll samples
# ==============================================================================

# ==============================================================================
# TIME-MATCHED APPROACH - Replaces daily mean approach
# ==============================================================================

cat("Processing discrete chlorophyll samples with time matching...\n")

# Step 1: Process bulk chlorophyll (keep collected time)
c_qc <- c %>% 
  select(date, collected, station = site_id, line_out_depth, filter_type, chla, chla_flag) %>% 
  filter(filter_type == "Bulk GF/F") %>% 
  filter(chla_flag %in% c("AV", "SVC", "ADL") | is.na(chla_flag)) %>%
  mutate(collected = lubridate::with_tz(collected, tzone = "America/Vancouver"))

# Prep bulk data
c_bulk <- c_qc %>% 
  mutate(pres = round(line_out_depth)) %>%
  mutate(pres = case_when(pres == 0 ~ 1, TRUE ~ as.numeric(pres))) %>%
  drop_na(chla) %>%
  select(date, collected, station, pres, chl_bulk = chla)

# Step 2: Process size-fractionated chlorophyll
c_sf <- c %>% 
  select(date, collected, station = site_id, line_out_depth, filter_type, chla, chla_flag) %>% 
  filter(!filter_type == "Bulk GF/F") %>% 
  filter(chla_flag %in% c("AV", "SVC", "ADL") | is.na(chla_flag))

# Calculate size-fractionated sum (keep collected time)
c_sf_sum <- c_sf %>% 
  filter(!is.na(chla)) %>%
  filter(chla > 0) %>% 
  # Standardize filter types (2um = 3um)
  mutate(filter_type_std = case_when(
    filter_type == "2um" ~ "3um",
    TRUE ~ as.character(filter_type)
  )) %>%
  # Average by depth and filter type (in case of duplicates at same time)
  group_by(date, collected, station, line_out_depth, filter_type_std) %>% 
  summarise(avg_chla = mean(chla), .groups = "drop") %>%
  # Check that we have all 3 filter types
  group_by(date, collected, station, line_out_depth) %>% 
  mutate(n_filt = n_distinct(filter_type_std)) %>% 
  ungroup() %>% 
  filter(n_filt == 3) %>%  # Must have all 3 filter types
  # Calculate sum
  group_by(date, collected, station, line_out_depth) %>% 
  summarise(chl_sf_sum = sum(avg_chla), .groups = "drop") %>%
  mutate(pres = round(line_out_depth)) %>%
  mutate(pres = case_when(pres == 0 ~ 1, TRUE ~ as.numeric(pres))) %>%
  select(date, collected, station, pres, chl_sf_sum)

# Step 3: Combine bulk and size-fractionated (prefer bulk, use SF when missing)
c_combined <- c_bulk %>%
  full_join(c_sf_sum, by = c("date", "collected", "station", "pres")) %>%
  mutate(chl_comb = case_when(
    !is.na(chl_bulk) ~ chl_bulk,
    is.na(chl_bulk) & !is.na(chl_sf_sum) ~ chl_sf_sum,
    TRUE ~ NA_real_
  )) %>%
  drop_na(chl_comb) %>%
  # Filter for minimum number of depths per cast
  group_by(date, collected, station) %>%
  mutate(n_depths = n()) %>%
  ungroup() %>%
  filter(n_depths >= 3)  # Need at least 3 depths for regression

cat("Discrete chlorophyll prepared with", nrow(c_combined), "samples\n")

# Step 4: Prepare CTD data (keep individual casts, not daily means)
f_casts <- f %>%
  select(castpk, hakai_id, date, station, time, ctdNum, pres, flu) %>%
  # For each cast, get a single representative time (use median pressure time)
  group_by(castpk, date, station) %>%
  mutate(cast_time = median(time)) %>%  # Representative time for the cast
  ungroup()

cat("CTD data prepared with", n_distinct(f_casts$castpk), "unique casts\n")

# Step 5: Match discrete samples to CTD casts by closest time
# First create a table of unique casts with their times
cast_times <- f_casts %>%
  distinct(castpk, date, station, ctdNum, cast_time)

# For each discrete sample, find the closest CTD cast
matches <- c_combined %>%
  left_join(cast_times, by = c("date", "station"), relationship = "many-to-many") %>%
  mutate(time_diff_hours = abs(as.numeric(difftime(collected, cast_time, units = "hours")))) %>%
  # Keep only matches within 6 hours
  filter(time_diff_hours <= 6) %>%
  # For each discrete sample, keep the closest cast
  group_by(date, station, collected, pres) %>%
  slice_min(time_diff_hours, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Matched", nrow(matches), "discrete samples to CTD casts\n")
cat("Time differences: median =", round(median(matches$time_diff_hours), 2), 
    "hours, max =", round(max(matches$time_diff_hours), 2), "hours\n")

# Step 6: Join matched discrete samples with CTD fluorescence data
f_fit <- matches %>%
  left_join(
    f_casts %>% select(castpk, date, station, pres, flu, time),
    by = c("castpk", "date", "station", "pres")
  ) %>%
  drop_na(flu, chl_comb) %>%
  rename(f_dm = flu) %>%  # Keep same name for compatibility with existing code
  mutate(cast_id = paste(date, station, castpk, sep = "-")) %>%
  select(cast_id, castpk, date, station, ctdNum, pres, time, collected, 
         time_diff_hours, f_dm, chl_comb, chl_bulk, chl_sf_sum)

# Also keep full CTD profiles for plotting (rename to f_dm for compatibility)
f_dm <- f_casts %>%
  rename(f_dm = flu) %>%
  mutate(cast_id = paste(date, station, castpk, sep = "-")) %>%
  select(cast_id, castpk, date, station, ctdNum, pres, time, f_dm)

cat("Final dataset:", nrow(f_fit), "matched CTD-discrete pairs\n")
cat("Covering", n_distinct(f_fit$cast_id), "unique casts\n")

# Check for casts with sufficient data for regression
cast_summary <- f_fit %>%
  group_by(cast_id, date, station, ctdNum) %>%
  summarise(
    n_points = n(),
    time_diff_range = paste(round(min(time_diff_hours), 2), "-", 
                            round(max(time_diff_hours), 2), "hrs"),
    .groups = "drop"
  )

cat("\nCasts by number of matched points:\n")
print(table(cast_summary$n_points))
# ==============================================================================
# STEP 4: Create dataset for linear fits
# ==============================================================================

cat("Preparing dataset for calibration fits...\n")

# Creating dataset to apply linear fits between chlorophyll and fluorometry
# f_fit <- f_dm %>%
#   left_join(c_join, by = c("date", "station", "pres", "id")) %>% 
#   drop_na(chl_comb.x, f_dm) %>%
#   select(-chl_comb.y, -chl_dm.y) %>%
#   rename(chl_comb = chl_comb.x, chl_dm = chl_dm.x) %>%
#   unite(id, c(date, station), sep = "-", remove = F)

cat("Calibration dataset prepared.\n")

# ==============================================================================
# STEP 5: Calculate cast-by-cast regression models
# ==============================================================================

cat("Calculating cast-by-cast regression models...\n")

# Fit both models (with and without surface) and select the best
model <- f_fit %>% 
  filter(!is.na(chl_comb) & !is.na(f_dm)) %>%
  group_by(cast_id) %>%  # Use cast_id directly
  summarise({
    # Check data availability
    data_no_surf <- filter(pick(everything()), pres != 1)  # Updated from cur_data()
    data_with_surf <- pick(everything())  # Updated from cur_data()
    
    nobs_no_surf <- nrow(data_no_surf)
    nobs_with_surf <- nrow(data_with_surf)
    
    # Initialize variables
    can_fit_no_surf <- nobs_no_surf >= 2
    can_fit_with_surf <- nobs_with_surf >= 2
    
    # Fit model WITHOUT surface (if possible)
    if (can_fit_no_surf) {
      model_no_surf <- lm(chl_comb ~ f_dm, data = data_no_surf)
      intercept_no_surf <- round(coef(model_no_surf)[1], 2)
      slope_no_surf <- round(coef(model_no_surf)[2], 2)
      r2_no_surf <- round(summary(model_no_surf)$adj.r.squared, 2)
      p.value_no_surf <- round(summary(model_no_surf)$coefficients[2, 4], 5)
    } else {
      intercept_no_surf <- slope_no_surf <- r2_no_surf <- p.value_no_surf <- NA
    }
    
    # Fit model WITH surface (if possible)
    if (can_fit_with_surf) {
      model_with_surf <- lm(chl_comb ~ f_dm, data = data_with_surf)
      intercept_with_surf <- round(coef(model_with_surf)[1], 2)
      slope_with_surf <- round(coef(model_with_surf)[2], 2)
      r2_with_surf <- round(summary(model_with_surf)$adj.r.squared, 2)
      p.value_with_surf <- round(summary(model_with_surf)$coefficients[2, 4], 5)
    } else {
      intercept_with_surf <- slope_with_surf <- r2_with_surf <- p.value_with_surf <- NA
    }
    
    # Select best model
    if (!can_fit_no_surf & !can_fit_with_surf) {
      # Neither model can be fit
      use_surface <- NA
      intercept <- slope <- r2 <- p.value <- nobs <- NA
    } else if (!can_fit_no_surf) {
      # Only surface model can be fit
      use_surface <- TRUE
      intercept <- intercept_with_surf
      slope <- slope_with_surf
      r2 <- r2_with_surf
      p.value <- p.value_with_surf
      nobs <- nobs_with_surf
    } else if (!can_fit_with_surf) {
      # Only no-surface model can be fit
      use_surface <- FALSE
      intercept <- intercept_no_surf
      slope <- slope_no_surf
      r2 <- r2_no_surf
      p.value <- p.value_no_surf
      nobs <- nobs_no_surf
    } else {
      # Both models can be fit - choose best based on p-value
      use_surface <- p.value_with_surf < p.value_no_surf
      intercept <- ifelse(use_surface, intercept_with_surf, intercept_no_surf)
      slope <- ifelse(use_surface, slope_with_surf, slope_no_surf)
      r2 <- ifelse(use_surface, r2_with_surf, r2_no_surf)
      p.value <- ifelse(use_surface, p.value_with_surf, p.value_no_surf)
      nobs <- ifelse(use_surface, nobs_with_surf, nobs_no_surf)
    }
    
    # Set slope to NA if p-value > 0.05 OR if slope is negative
    slope_best <- ifelse(is.na(p.value) || p.value > 0.05 || slope < 0, NA, slope)
    
    # Determine if surface is quenched (available but not used)
    surface_available <- any(data_with_surf$pres == 1)
    quenched <- surface_available && !isTRUE(use_surface)
    
    data.frame(
      intercept = intercept,
      slope = slope_best,
      r2 = r2,
      p.value = p.value,
      nobs = nobs,
      use_surface = use_surface,
      quenched = quenched
    )
  }, .groups = 'drop') %>%
  # Don't separate cast_id - extract date and station directly from f_fit
  left_join(
    f_fit %>% distinct(cast_id, date, station, castpk, ctdNum),
    by = "cast_id"
  ) %>%
  arrange(date, station, castpk)

cat("Regression models calculated for", nrow(model), "casts\n")

# ==============================================================================
# STEP 6: Calculate sensor-specific statistics and flag outliers
# ==============================================================================

cat("Calculating sensor-specific statistics...\n")

# Actually, since we already joined ctdNum in the model calculation above,
# this should be even simpler:
model_with_sensor <- model  # ctdNum is already there!

# Calculate sensor medians for reference lines
sensor_medians <- model_with_sensor %>%
  filter(!is.na(slope)) %>%
  group_by(ctdNum) %>%
  summarise(
    median_slope = median(slope, na.rm = TRUE),
    median_intercept = median(intercept, na.rm = TRUE)
  )

# Flag casts with unusual slopes or poor fits - BY SENSOR
model_flagged <- model_with_sensor %>%
  group_by(ctdNum) %>%
  mutate(
    # Flag slopes far from median WITHIN each sensor
    median_slope = median(slope, na.rm = TRUE),
    mad_slope = mad(slope, na.rm = TRUE),
    slope_z = abs(slope - median_slope) / mad_slope,
    flag_slope = slope_z > 3,
    
    # Flag poor fits
    flag_r2 = r2 < 0.7,
    flag_pvalue = p.value > 0.01,
    
    # Flag unusual intercepts
    median_intercept = median(intercept, na.rm = TRUE),
    mad_intercept = mad(intercept, na.rm = TRUE),
    intercept_z = abs(intercept - median_intercept) / mad_intercept,
    flag_intercept = intercept_z > 3,
    
    # Auto-flag from algorithm
    auto_flag = flag_slope | flag_r2 | flag_intercept
  ) %>%
  ungroup()

# Create flagged summary
flagged_summary <- model_flagged %>%
  filter(auto_flag) %>%
  select(ctdNum, cast_id, date, station, slope, median_slope, slope_z,
         intercept, r2, p.value, nobs, 
         flag_slope, flag_r2, flag_intercept) %>%
  arrange(ctdNum, desc(slope_z))

cat("Auto-flagging complete. Found", nrow(flagged_summary), "auto-flagged casts.\n")

# ==============================================================================
# STEP 7: Initialize QC tracking dataframe
# ==============================================================================

# Check if QC file already exists
qc_file <- here("figures", "cast_qc_flags.csv")

# Always rebuild cast_qc to match current cast_ids
cast_qc <- model_flagged %>%
  select(cast_id, ctdNum, date, station, auto_flag) %>%
  mutate(
    user_flag = FALSE,
    qc_comment = "",
    qc_date = as.character(NA)
  )

# If old QC file exists, try to merge in previous user flags/comments
if (file.exists(qc_file)) {
  old_qc <- read.csv(qc_file, stringsAsFactors = FALSE)
  cat("Found existing QC file. Attempting to merge previous flags...\n")
  
  # Try to match by date-station if cast_id format changed
  if ("cast_id" %in% names(old_qc)) {
    # Extract date-station from old cast_ids for matching
    old_qc <- old_qc %>%
      mutate(
        date_station = sub("-[^-]+$", "", cast_id)  # Remove last segment (castpk)
      ) %>%
      filter(user_flag == TRUE | nchar(qc_comment) > 0)  # Only keep user-modified entries
    
    cast_qc <- cast_qc %>%
      mutate(date_station = paste(date, station, sep = "-")) %>%
      left_join(
        old_qc %>% select(date_station, old_user_flag = user_flag, 
                          old_qc_comment = qc_comment, old_qc_date = qc_date),
        by = "date_station"
      ) %>%
      mutate(
        user_flag = ifelse(!is.na(old_user_flag), old_user_flag, user_flag),
        qc_comment = ifelse(!is.na(old_qc_comment), old_qc_comment, qc_comment),
        qc_date = ifelse(!is.na(old_qc_date), old_qc_date, qc_date)
      ) %>%
      select(-date_station, -old_user_flag, -old_qc_comment, -old_qc_date)
    
    cat("Merged", sum(cast_qc$user_flag, na.rm = TRUE), "user flags from previous QC file\n")
  }
} else {
  cat("Initialized new QC tracking dataframe\n")
}

cat("QC tracking ready with", nrow(cast_qc), "casts\n")

# ==============================================================================
# STEP 8: Plotting function
# ==============================================================================

plot_cast_detailed <- function(cast_id) {
  # Parse date, station, and castpk from cast_id
  parts <- strsplit(cast_id, "-")[[1]]
  cast_date <- as.Date(paste(parts[1:3], collapse = "-"))
  cast_station <- parts[4]
  cast_castpk <- parts[5]
  
  # Filter by cast_id directly
  profile_data <- f_dm %>%
    filter(cast_id == !!cast_id)
  
  if (nrow(profile_data) == 0) {
    stop(paste("No profile data found for cast:", cast_id))
  }
  
  # Get discrete match data for scatter plot
  cast_data <- f_fit %>%
    filter(cast_id == !!cast_id)
  
  # Get sensor info
  sensor <- profile_data$ctdNum[1]
  sensor_med <- sensor_medians %>% filter(ctdNum == sensor)
  
  # Get cast info - EXTRACT SINGLE ROW
  cast_info <- model_flagged %>% 
    filter(cast_id == !!cast_id) %>%
    slice(1)  # Ensure single row
  
  # Get QC status - EXTRACT SINGLE ROW
  qc_status <- cast_qc %>% 
    filter(cast_id == !!cast_id) %>%
    slice(1)  # Ensure single row
  
  # Check if we have the info we need
  if (nrow(cast_info) == 0) {
    stop(paste("No model info found for cast:", cast_id))
  }
  if (nrow(qc_status) == 0) {
    stop(paste("No QC status found for cast:", cast_id))
  }
  
  # Fit model (potentially excluding surface)
  if (cast_info$use_surface[1]) {  # Extract single value
    fit_data <- cast_data %>% filter(!is.na(chl_comb), !is.na(f_dm))
  } else {
    fit_data <- cast_data %>% filter(!is.na(chl_comb), !is.na(f_dm), pres != 1)
  }
  
  if (nrow(fit_data) < 2) {
    stop(paste("Not enough data points for regression in cast:", cast_id))
  }
  
  cast_model <- lm(chl_comb ~ f_dm, data = fit_data)
  cast_slope <- coef(cast_model)[2]
  cast_intercept <- coef(cast_model)[1]
  r2 <- summary(cast_model)$r.squared
  p_val <- summary(cast_model)$coefficients[2, 4]
  
  # Equation text
  eq_label <- sprintf(
    "Cast: y = %.2fx + %.2f\nMedian: y = %.2fx + %.2f\nR² = %.2f, p = %.4f",
    cast_slope, cast_intercept,
    sensor_med$median_slope[1], sensor_med$median_intercept[1],
    r2, p_val
  )
  
  # Flag status text - EXTRACT SINGLE VALUES
  flag_text <- paste(
    ifelse(qc_status$auto_flag[1], "AUTO-FLAGGED", ""),
    ifelse(qc_status$user_flag[1], "USER-FLAGGED", ""),
    sep = " "
  )
  if (flag_text == " ") flag_text <- "NOT FLAGGED"
  
  profile_plot <- profile_data %>%
    filter(pres < 50) %>%
    ggplot(aes(y = pres)) +
    geom_line(aes(x = f_dm, color = "CTD Fluorescence"),
              linewidth = 1, orientation = "y") +
    geom_point(data = cast_data %>% filter(pres < 50),
               aes(x = chl_comb, color = "Discrete Samples"), 
               size = 3, shape = 16) +
    scale_y_reverse(name = "Pressure (dbar)") +
    scale_x_continuous(name = "Chlorophyll (µg/L)") +
    scale_color_manual(
      name = "",
      values = c("CTD Fluorescence" = "#2E86AB", 
                 "Discrete Samples" = "#A23B72"),
      breaks = c("CTD Fluorescence", "Discrete Samples")
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    ) +
    labs(title = "A) Fluorescence Profile")
  
  # Calculate axis limits for scatter
  axis_min <- min(c(cast_data$f_dm, cast_data$chl_comb), na.rm = TRUE)
  axis_max <- max(c(cast_data$f_dm, cast_data$chl_comb), na.rm = TRUE)
  buffer <- (axis_max - axis_min) * 0.05
  axis_limits <- c(max(0, axis_min - buffer), axis_max + buffer)
  
  # SCATTER PLOT - only discrete matches
  scatter_plot <- cast_data %>%
    filter(!is.na(chl_comb), !is.na(f_dm)) %>%
    ggplot(aes(x = f_dm, y = chl_comb, color = as.factor(pres))) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_abline(slope = sensor_med$median_slope[1], 
                intercept = sensor_med$median_intercept[1],
                color = "blue", linewidth = 1.2, linetype = "solid") +
    geom_smooth(data = fit_data,
                aes(x = f_dm, y = chl_comb),
                method = "lm", se = TRUE,
                color = "red", fill = "red",
                linewidth = 1, alpha = 0.2,
                inherit.aes = FALSE) +
    geom_point(size = 3) +
    scale_color_brewer(
      name = "Pressure\n(dbar)",
      palette = "Set2"
    ) +
    scale_x_continuous(name = "CTD Fluorescence (µg/L)", 
                       limits = axis_limits) +
    scale_y_continuous(name = "Discrete Chlorophyll (µg/L)", 
                       limits = axis_limits) +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    ) +
    labs(
      title = "B) CTD vs Discrete Chlorophyll",
      subtitle = paste0(
        "Red = cast fit | Blue = sensor median | Dashed = 1:1",
        ifelse(!cast_info$use_surface[1], " | Surface excluded", "")  # Extract single value
      )
    ) +
    annotate("text", 
             x = axis_limits[1], 
             y = axis_limits[2], 
             label = eq_label,
             hjust = -0.05, 
             vjust = 1.1,
             size = 3.5,
             fontface = "italic")
  
  # COMBINE PLOTS
  combined_plot <- profile_plot + scatter_plot +
    plot_annotation(
      title = paste0("Cast: ", cast_id, " | Sensor: ", sensor, " | ", flag_text),
      subtitle = if(nchar(qc_status$qc_comment[1]) > 0) {  # Extract single value
        paste0("Comment: ", qc_status$qc_comment[1])
      } else {
        NULL
      },
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, color = "red")
      )
    )
  
  return(combined_plot)
}

# ==============================================================================
# STEP 9: Interactive browsing function - organized by sensor
# ==============================================================================

browse_casts <- function(flagged_only = FALSE) {
  # Create subfolder for QC flagged casts
  qc_figures_dir <- here("figures", "qc_flagged_casts")
  if (!dir.exists(qc_figures_dir)) {
    dir.create(qc_figures_dir, recursive = TRUE)
    cat("Created directory:", qc_figures_dir, "\n")
  }
  
  # Get cast list organized by sensor and date
  if (flagged_only) {
    cast_df <- model_flagged %>% 
      filter(auto_flag) %>% 
      arrange(ctdNum, date) %>%
      select(cast_id, ctdNum, date, station)
  } else {
    cast_df <- model_flagged %>% 
      filter(!is.na(slope)) %>%
      arrange(ctdNum, date) %>%
      select(cast_id, ctdNum, date, station)
  }
  
  if (nrow(cast_df) == 0) {
    cat("No casts to browse!\n")
    return(invisible())
  }
  
  cast_list <- cast_df$cast_id
  
  # Get sensor info for display
  sensor_summary <- cast_df %>%
    group_by(ctdNum) %>%
    summarise(n_casts = n())
  
  cat(paste0("\n=== Browsing ", nrow(cast_df), " casts organized by sensor ===\n"))
  cat("Sensors and cast counts:\n")
  print(sensor_summary)
  cat("\n")
  
  i <- 1
  current_sensor <- NULL
  
  repeat {
    # Bounds checking
    if (i < 1) i <- 1
    if (i > length(cast_list)) {
      cat("\nEnd of cast list.\n")
      break
    }
    
    cast_id <- cast_list[i]
    cast_sensor <- cast_df$ctdNum[i]
    
    # Announce when switching sensors
    if (is.null(current_sensor) || current_sensor != cast_sensor) {
      sensor_casts <- sum(cast_df$ctdNum == cast_sensor)
      sensor_position <- which(cast_df$ctdNum == cast_sensor)[1]
      cat(paste0("\n", strrep("=", 70), "\n"))
      cat(paste0(">>> NOW BROWSING SENSOR ", cast_sensor, " (", sensor_casts, " casts) <<<\n"))
      cat(paste0(strrep("=", 70), "\n\n"))
      current_sensor <- cast_sensor
    }
    
    # Plot
    tryCatch({
      p <- plot_cast_detailed(cast_id)
      print(p)
      
      # Show current QC status - EXTRACT SINGLE ROW
      qc_status <- cast_qc %>% 
        filter(cast_id == !!cast_id) %>%
        slice(1)
      
      cat(paste0("\n=== Cast ", i, " of ", length(cast_list), ": ", cast_id, " ===\n"))
      cat(paste0("Sensor: ", cast_sensor, "\n"))
      cat(paste0("Auto-flag: ", ifelse(is.na(qc_status$auto_flag[1]), FALSE, qc_status$auto_flag[1]), "\n"))
      cat(paste0("User-flag: ", ifelse(is.na(qc_status$user_flag[1]), FALSE, qc_status$user_flag[1]), "\n"))
      if (!is.na(qc_status$qc_comment[1]) && nchar(qc_status$qc_comment[1]) > 0) {
        cat(paste0("Comment: ", qc_status$qc_comment[1], "\n"))
      }
    }, error = function(e) {
      cat("Error plotting cast:", cast_id, "\n")
      cat(e$message, "\n")
    })
    
    # Get user input
    cat("\nCommands:\n")
    cat("  [n]ext | [p]revious | [f]lag/unflag | [c]omment | [q]uit\n")
    cat("  Or enter cast number to jump\n")
    cat("Choice: ")
    input <- readline()
    
    if (input == "n" || input == "") {
      i <- i + 1
    } else if (input == "p") {
      i <- i - 1
      # Reset current_sensor when going backward to re-announce if needed
      if (i > 0) {
        prev_sensor <- cast_df$ctdNum[i]
        if (prev_sensor != current_sensor) {
          current_sensor <- NULL
        }
      }
    } else if (input == "f") {
      # Toggle flag - get index of matching row
      row_idx <- which(cast_qc$cast_id == cast_id)[1]
      
      if (is.na(row_idx)) {
        cat("Error: cast not found in QC tracking\n")
        i <- i  # Stay on same cast
        next
      }
      
      current_flag <- cast_qc$user_flag[row_idx]
      # Handle NA as FALSE
      if (is.na(current_flag)) current_flag <- FALSE
      
      cast_qc$user_flag[row_idx] <<- !current_flag
      cast_qc$qc_date[row_idx] <<- as.character(Sys.Date())
      cat(paste0("User flag set to: ", !current_flag, "\n"))
      
      # Auto-save if flagged or has comment - HANDLE NAs
      qc_row <- cast_qc[row_idx, ]
      user_flagged <- isTRUE(qc_row$user_flag[1])
      has_comment <- !is.na(qc_row$qc_comment[1]) && nchar(qc_row$qc_comment[1]) > 0
      
      if (user_flagged || has_comment) {
        filename <- paste0("cast_", gsub("-", "_", cast_id), ".png")
        ggsave(file.path(qc_figures_dir, filename),
               plot = p, width = 12, height = 6, dpi = 300, bg = "white")
        cat(paste0("→ Auto-saved to: ", file.path(qc_figures_dir, filename), "\n"))
      }
      
      # Re-plot to show updated flag
      i <- i  # Stay on same cast
    } else if (input == "c") {
      cat("Enter comment (or press Enter to clear): ")
      comment <- readline()
      
      # Get index of matching row
      row_idx <- which(cast_qc$cast_id == cast_id)[1]
      
      if (is.na(row_idx)) {
        cat("Error: cast not found in QC tracking\n")
        i <- i  # Stay on same cast
        next
      }
      
      cast_qc$qc_comment[row_idx] <<- comment
      cast_qc$qc_date[row_idx] <<- as.character(Sys.Date())
      cat("Comment saved.\n")
      
      # Auto-save if flagged or has comment - HANDLE NAs
      qc_row <- cast_qc[row_idx, ]
      user_flagged <- isTRUE(qc_row$user_flag[1])
      has_comment <- !is.na(qc_row$qc_comment[1]) && nchar(qc_row$qc_comment[1]) > 0
      
      if (user_flagged || has_comment) {
        filename <- paste0("cast_", gsub("-", "_", cast_id), ".png")
        ggsave(file.path(qc_figures_dir, filename),
               plot = p, width = 12, height = 6, dpi = 300, bg = "white")
        cat(paste0("→ Auto-saved to: ", file.path(qc_figures_dir, filename), "\n"))
      }
      
      # Re-plot to show updated comment
      i <- i  # Stay on same cast
    } else if (input == "q") {
      # Save QC data before quitting
      write.csv(cast_qc, qc_file, row.names = FALSE)
      cat(paste0("\nQC flags saved to: ", qc_file, "\n"))
      break
    } else if (grepl("^[0-9]+$", input)) {
      i <- as.numeric(input)
      # Reset current_sensor to re-announce when jumping
      current_sensor <- NULL
    } else {
      cat("Invalid input. Try again.\n")
    }
  }
  
  # Final save
  write.csv(cast_qc, qc_file, row.names = FALSE)
  cat(paste0("\nFinal QC flags saved to: ", qc_file, "\n"))
  
  # Summary by sensor
  summary_by_sensor <- cast_qc %>%
    left_join(cast_df %>% select(cast_id, ctdNum), by = "cast_id") %>%
    group_by(ctdNum) %>%
    summarise(
      total_casts = n(),
      auto_flagged = sum(auto_flag, na.rm = TRUE),
      user_flagged = sum(user_flag, na.rm = TRUE),
      with_comments = sum(nchar(qc_comment) > 0, na.rm = TRUE)
    )
  
  cat("\n=== QC Summary by Sensor ===\n")
  print(summary_by_sensor)
  
  # Overall summary
  summary_stats <- cast_qc %>%
    summarise(
      total_casts = n(),
      auto_flagged = sum(auto_flag, na.rm = TRUE),
      user_flagged = sum(user_flag, na.rm = TRUE),
      with_comments = sum(nchar(qc_comment) > 0, na.rm = TRUE)
    )
  
  cat("\n=== Overall QC Summary ===\n")
  print(summary_stats)
  cat(paste0("\nFlagged cast figures saved to: ", qc_figures_dir, "\n"))
}

# ==============================================================================
# STEP 10: Export final flagged list
# ==============================================================================

export_qc_summary <- function() {
  qc_summary <- cast_qc %>%
    left_join(model_flagged %>% select(cast_id, slope, r2, p.value, slope_z), 
              by = "cast_id") %>%
    filter(auto_flag | user_flag) %>%
    arrange(ctdNum, date)
  
  output_file <- here("figures", "cast_qc_summary.csv")
  write.csv(qc_summary, output_file, row.names = FALSE)
  cat(paste0("QC summary exported to: ", output_file, "\n"))
  
  return(qc_summary)
}

# ==============================================================================
# USAGE INSTRUCTIONS
# ==============================================================================
cat("\n=== CTD Fluorescence QC Tool Ready ===\n")
cat("\nTo start:\n")
cat("  browse_casts()              # Browse all casts\n")
cat("  browse_casts(flagged_only = TRUE)  # Browse only auto-flagged casts\n")
cat("\nTo export final summary:\n")
cat("  export_qc_summary()\n")
cat("\nQC data will be saved to: ", qc_file, "\n\n")




