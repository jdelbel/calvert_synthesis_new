library(dplyr)
library(ggplot2)
library(lubridate)
library(zoo)
library(tidyverse)
library(here)

# ============================================================================
# FLUOROMETRY QUALITY CONTROL PIPELINE
# ============================================================================

# Data preparation function with flexible column names and gap handling
prepare_fluorometry_data <- function(chl_data, 
                                     years = NULL,
                                     bad_periods = NULL,
                                     fluorescence_col = c("Chlorophyll", "Fluorometer")) {
  
  cat("=== DATA PREPARATION ===\n")
  
  # Filter years if specified
  if(!is.null(years)) {
    chl_data <- chl_data %>% filter(year %in% years)
    cat("Filtered to years:", paste(years, collapse = ", "), "\n")
  }
  
  # Check which fluorescence columns are available
  available_fl_cols <- intersect(fluorescence_col, names(chl_data))
  
  if(length(available_fl_cols) == 0) {
    stop("No fluorescence column found. Available columns: ", paste(names(chl_data), collapse = ", "))
  }
  
  cat("Available fluorescence columns:", paste(available_fl_cols, collapse = ", "), "\n")
  
  # Step 1: Create datetime columns first
  data_prep <- chl_data %>%
    mutate(
      measurementTime_2 = lubridate::mdy_hm(measurementTime),
      date = as.Date(lubridate::date(measurementTime_2)),
      month = lubridate::month(measurementTime_2),
      hour = lubridate::hour(measurementTime_2),
      year = lubridate::year(measurementTime_2)
    )
  
  # Step 2: Create unified fluorescence column
  if(length(available_fl_cols) == 1) {
    # Only one column available
    data_prep <- data_prep %>%
      mutate(fl = !!sym(available_fl_cols[1]))
    cat("Using single fluorescence column:", available_fl_cols[1], "\n")
    
  } else {
    # Multiple columns available - combine them
    cat("Combining multiple fluorescence columns with coalesce():\n")
    
    data_prep <- data_prep %>%
      mutate(fl = coalesce(!!!syms(available_fl_cols)))
    
    # Report which columns had data
    for(col in available_fl_cols) {
      n_values <- sum(!is.na(chl_data[[col]]), na.rm = TRUE)
      if(n_values > 0) {
        cat("  ", col, ":", n_values, "values\n")
      }
    }
  }
  
  # Step 3: Final cleanup and arrange
  data_prep <- data_prep %>%
    select(measurementTime, fl, measurementTime_2, date, month, hour, year) %>%
    arrange(measurementTime_2)
  
  # Handle bad data periods and create groups
  if(!is.null(bad_periods)) {
    cat("Processing", length(bad_periods), "bad data periods\n")
    
    # Debug: Check date column class
    cat("Date column class:", class(data_prep$date), "\n")
    cat("Sample dates:", head(data_prep$date), "\n")
    
    # Ensure date is actually Date class
    data_prep <- data_prep %>%
      mutate(date = as.Date(date))
    
    cat("Date column class after conversion:", class(data_prep$date), "\n")
    
    # Initialize all data as group 1
    data_prep$group <- 1L  # Use integer
    current_group <- 1L
    
    # Process each bad period
    for(i in seq_along(bad_periods)) {
      period <- bad_periods[[i]]
      
      # Convert dates to Date class for consistent comparison
      period_start <- as.Date(period$start)
      period_end <- as.Date(period$end)
      
      cat(sprintf("  Bad period %d: %s to %s (%s)\n", 
                  i, period_start, period_end, 
                  ifelse(is.null(period$description), "No description", period$description)))
      
      cat("  Period start class:", class(period_start), "\n")
      cat("  Period end class:", class(period_end), "\n")
      
      # Remove data within bad periods first
      n_before <- nrow(data_prep)
      data_prep <- data_prep %>%
        filter(date < period_start | date > period_end)
      n_after <- nrow(data_prep)
      
      cat("  Removed", n_before - n_after, "points in bad period\n")
      
      # Increment group number after each bad period
      current_group <- current_group + 1L
      
      # Create logical vector for group assignment
      if(nrow(data_prep) > 0) {
        after_period <- data_prep$date > period_end
        if(any(after_period, na.rm = TRUE)) {
          data_prep$group[after_period] <- current_group
        }
      }
    }
    
    # Summary of groups
    group_summary <- data_prep %>%
      group_by(group, year) %>%
      summarise(
        start_date = min(date),
        end_date = max(date),
        n_points = n(),
        .groups = "drop"
      )
    
    cat("\nData groups created:\n")
    print(group_summary)
    
  } else {
    # If no bad periods specified, create simple yearly groups
    data_prep$group <- as.numeric(factor(data_prep$year))
    cat("No bad periods specified. Created yearly groups.\n")
  }
  
  cat("Final dataset:", nrow(data_prep), "points in", max(data_prep$group), "groups\n\n")
  
  return(data_prep)
}

# Iterative cleaning function
clean_fluorometry_iterative <- function(data, max_iter = 3, mad_threshold = 3) {
  cat("Starting iterative cleaning...\n")
  
  for(i in 1:max_iter) {
    n_before <- nrow(data)
    
    data <- data %>%
      group_by(group) %>%
      mutate(
        gf = c(NA, diff(fl)/head(fl, -1)),
        gf_mad = mad(gf, na.rm = TRUE),
        gf_median = median(gf, na.rm = TRUE),
        outlier = abs(gf - gf_median) > mad_threshold * gf_mad | 
          fl < 0 | 
          fl > 50
      ) %>%
      filter(!outlier) %>%
      select(-gf_mad, -gf_median, -outlier) %>%
      ungroup()
    
    n_after <- nrow(data)
    removed <- n_before - n_after
    
    cat(sprintf("  Iteration %d: removed %d points (%.2f%%)\n", 
                i, removed, (removed/n_before)*100))
    
    if(removed < 0.001 * n_before) break
  }
  return(data)
}

# Complete QC Pipeline Function
fluorometry_qc_pipeline <- function(raw_data) {
  
  cat("=== FLUOROMETRY QC PIPELINE ===\n")
  cat("Starting with", nrow(raw_data), "raw data points\n\n")
  
  # Store data at each step for plotting
  qc_steps <- list()
  qc_steps$raw <- raw_data
  
  # ---- STEP 1: Iterative spike removal ----
  cat("STEP 1: Iterative spike removal\n")
  step1 <- clean_fluorometry_iterative(raw_data)
  cat("  Retained:", nrow(step1), "points\n")
  qc_steps$step1_spikes <- step1
  
  # ---- STEP 2: Nighttime filtering ----
  cat("\nSTEP 2: Nighttime filtering (removing quenched data)\n")
  step2 <- step1 %>%
    group_by(date) %>%
    filter(hour <= 6 | hour >= 20) %>%
    ungroup() %>%
    mutate(date_corr = case_when(hour >= 0 & hour <= 6 ~ (date-1),
                                 TRUE ~ as.Date(date)))
  
  cat("  Retained:", nrow(step2), "points\n")
  qc_steps$step2_nighttime <- step2
  
  # ---- STEP 3: Z-score filtering ----
  cat("\nSTEP 3: Z-score outlier removal (|z| > 3)\n")
  step3 <- step2 %>%
    group_by(month, group) %>%
    mutate(sd_all = sd(fl),
           mean_all = mean(fl),
           sd_mean_all = (fl - mean_all)/sd_all,
           cv = (sd_all/mean_all)*100) %>%
    filter(sd_mean_all < 3 & sd_mean_all > -3) %>%
    ungroup()
  
  cat("  Retained:", nrow(step3), "points\n")
  qc_steps$step3_zscore <- step3
  
  # ---- STEP 4: Daily medians ----
  cat("\nSTEP 4: Calculate daily medians\n")
  final <- step3 %>%
    group_by(date_corr) %>%
    mutate(fl_med_day = median(fl), 
           fl_med_3 = zoo::rollapply(fl, 3, median, align = 'right', fill = NA)) %>%
    ungroup()
  
  qc_steps$final <- final
  
  cat("\n=== PIPELINE COMPLETE ===\n")
  cat("Final dataset:", nrow(final), "points\n")
  
  return(qc_steps)
}

# Plotting function
plot_qc_progression <- function(qc_steps) {
  
  # Prepare data for plotting
  plot_data <- list()
  
  # Raw data
  plot_data$`Raw Data` <- qc_steps$raw %>%
    select(measurementTime_2, fl)
  
  # Step 1: Spike removal
  plot_data$`After Spike Removal` <- qc_steps$step1_spikes %>%
    select(measurementTime_2, fl)
  
  # Step 2: Nighttime only
  plot_data$`Nighttime Only` <- qc_steps$step2_nighttime %>%
    select(measurementTime_2, fl)
  
  # Step 3: Z-score filtered
  plot_data$`Z-score Filtered` <- qc_steps$step3_zscore %>%
    select(measurementTime_2, fl)
  
  # Step 4: Daily medians
  plot_data$`Daily Medians` <- qc_steps$final %>%
    select(measurementTime_2, fl_med_day) %>%
    rename(fl = fl_med_day)
  
  # Combine all data
  combined_data <- bind_rows(plot_data, .id = "step")
  
  # Create factor levels for proper ordering
  step_order <- names(plot_data)
  combined_data$step <- factor(combined_data$step, levels = step_order)
  
  # Create the plot
  p <- ggplot(combined_data, aes(x = measurementTime_2, y = fl)) +
    geom_line(size = 0.3, color = "black") +
    facet_wrap(~step, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(
      strip.text = element_text(hjust = 0, face = "bold", size = 10),
      panel.spacing = unit(0.1, "lines"),
      axis.text.x = element_text(angle = 0),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    labs(
      x = "Date",
      y = "Fluorescence",
      title = "Fluorometry QC Pipeline Progression"
    ) +
    scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y")
  
  return(p)
}

# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# Example 1: Single year (2024) with bad data period
# bad_periods_2024 <- list(
#   list(start = "2024-03-20", end = "2024-03-25", description = "Instrument failure")
# )
# 
# # Prepare data
# chl_22 <- prepare_fluorometry_data(
#   chl_data = chl, 
#   years = 2024,
#   bad_periods = bad_periods_2024,
#   fluorescence_col = c("Chlorophyll", "Fluorometer")
# )
# 
# # Run QC pipeline
# qc_results <- fluorometry_qc_pipeline(chl_22)
# progression_plot <- plot_qc_progression(qc_results)
# print(progression_plot)

# Example 2: Multiple years with various bad periods
# bad_periods_multi <- list(
#   list(start = "2022-07-15", end = "2022-07-20", description = "Buoy turnaround"),
#   list(start = "2023-05-01", end = "2023-05-10", description = "Maintenance"),
#   list(start = "2024-03-20", end = "2024-03-25", description = "Instrument failure")
# )
# 
# # Prepare multi-year data
# chl_multi <- prepare_fluorometry_data(
#   chl_data = chl,
#   years = c(2022, 2023, 2024),
#   bad_periods = bad_periods_multi
# )
# 
# # Run QC pipeline
# qc_results_multi <- fluorometry_qc_pipeline(chl_multi)

# Example 3: No bad periods (simple yearly grouping)
# chl_simple <- prepare_fluorometry_data(
#   chl_data = chl,
#   years = c(2023, 2024),
#   bad_periods = NULL
# )








#Downloading hourly data from the KC10 buoy fl
chl <- read_csv(here("files", "2025-08-20.1hourSamples.all.csv"))

# 2024-12-05.1hourSamples.all
# 2025-08-20.1hourSamples.all

# Define your 2024 bad period
bad_periods_2024 <- list(
  list(start = "2023-05-06", end = "2023-12-31", description = "Instrument failure"),
  list(start = "2024-03-20", end = "2024-03-25", description = "Instrument failure"),
  list(start = "2024-12-24", end = "2024-12-31", description = "Instrument failure"),
  list(start = "2025-04-28", end = "2025-05-11", description = "Turn Around"),
  list(start = "2025-07-10", end = "2025-07-14", description = "Instrument failure")
)

# # Prepare data
chl_multi <- prepare_fluorometry_data(
  chl_data = chl,
  years = c(2022, 2023, 2024, 2025),
  bad_periods = bad_periods_2024
)

  # # Run the pipeline
qc_results <- fluorometry_qc_pipeline(chl_multi)
  # 
  # # Create progression plot
progression_plot <- plot_qc_progression(qc_results)
print(progression_plot)
  # 
# # Access final cleaned dataset
final_data <- qc_results$final
  
final_data %>% 
  select(date_corr, fl_med_day, group) %>% 
  distinct() %>% 
  filter(date_corr > "2022-01-01") %>% 
  ggplot(aes(x = date_corr, y = fl_med_day, color = as.factor(group))) +
  geom_line()

write.csv(final_data, here("outputs", "qc_buoy_2025-08-20.csv."))
