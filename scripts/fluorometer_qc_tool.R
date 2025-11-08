

library(tidyverse)
library(here)

# Function to run the QC flagging
flag_fluorescence_profiles <- function(prof) {
  
  # Get unique cast identifiers with date info
  cast_summary <- prof %>%
    filter(pres <= 100) %>%
    group_by(castpk) %>%
    summarise(
      date = first(date),
      station = first(station),
      max_flu = max(flu, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(date, castpk)
  
  # Initialize results dataframe
  results <- data.frame(
    castpk = character(),
    date = as.Date(character()),
    flag = character(),
    comment = character(),
    stringsAsFactors = FALSE
  )
  
  # Check if previous results exist
  date_str <- format(Sys.Date(), "%Y%m%d")
  output_file <- here("outputs", paste0("fluorescence_QC_flags_", date_str, ".csv"))
  
  if (file.exists(output_file)) {
    cat("\nPrevious flagging file found. Loading...\n")
    results <- read_csv(output_file, show_col_types = FALSE)
  }
  
  # Get already flagged casts
  flagged_casts <- results$castpk
  
  # Filter to unflagged casts
  remaining_casts <- cast_summary %>%
    filter(!castpk %in% flagged_casts)
  
  cat("\n=== FLUORESCENCE PROFILE QC TOOL ===\n")
  cat(sprintf("Total casts: %d\n", nrow(cast_summary)))
  cat(sprintf("Already flagged: %d\n", length(flagged_casts)))
  cat(sprintf("Remaining: %d\n\n", nrow(remaining_casts)))
  
  if (nrow(remaining_casts) == 0) {
    cat("All casts have been flagged!\n")
    return(results)
  }
  
  cat("KEYBOARD SHORTCUTS:\n")
  cat("  1 or A = AV (Accepted Value)\n")
  cat("  2 or C = SVC (Suspicious Value - Caution)\n")
  cat("  3 or D = SVD (Suspicious Value - Do Not Use)\n")
  cat("  4 or S = SAT (Saturated)\n")
  cat("  B = Go back to previous cast\n")
  cat("  Q = Quit and save\n\n")
  
  i <- 1
  previous_flags <- list()
  
  while (i <= nrow(remaining_casts)) {
    
    current_cast <- remaining_casts$castpk[i]
    current_date <- remaining_casts$date[i]
    current_station <- remaining_casts$station[i]
    
    # Get profile data
    profile_data <- prof %>%
      filter(castpk == current_cast, pres <= 100) %>%
      arrange(pres)
    
    # Plot the profile
    plot(profile_data$flu, -profile_data$pres,
         type = "l", lwd = 2, col = "darkblue",
         xlab = "Fluorescence (ug/L)",
         ylab = "Pressure (dbar)",
         main = sprintf("Cast: %s | Date: %s | Station: %s\n[%d of %d]",
                        current_cast, current_date, current_station,
                        i, nrow(remaining_casts)),
         cex.main = 0.9)
    grid()
    points(profile_data$flu, -profile_data$pres, pch = 20, col = "darkblue")
    
    # Get user input
    cat(sprintf("\nCast %s (%d/%d) - Enter flag (1/A, 2/C, 3/D, 4/S, B=back, Q=quit): ",
                current_cast, i, nrow(remaining_casts)))
    
    user_input <- toupper(trimws(readline()))
    
    if (user_input == "Q") {
      cat("\nQuitting and saving...\n")
      break
    } else if (user_input == "B" && i > 1) {
      # Go back - remove last entry from results
      if (length(previous_flags) > 0) {
        last_flag <- previous_flags[[length(previous_flags)]]
        results <- results %>% filter(castpk != last_flag$castpk)
        previous_flags <- previous_flags[-length(previous_flags)]
      }
      i <- i - 1
      next
    } else if (user_input %in% c("1", "A", "2", "C", "3", "D", "4", "S")) {
      
      # Map input to flag
      flag <- case_when(
        user_input %in% c("1", "A") ~ "AV",
        user_input %in% c("2", "C") ~ "SVC",
        user_input %in% c("3", "D") ~ "SVD",
        user_input %in% c("4", "S") ~ "SAT",
        TRUE ~ NA_character_
      )
      
      # Get comment if not AV
      comment <- ""
      if (flag != "AV") {
        cat("Enter comment (or press Enter to skip): ")
        comment <- readline()
      }
      
      # Add to results
      new_row <- data.frame(
        castpk = current_cast,
        date = current_date,
        flag = flag,
        comment = comment,
        stringsAsFactors = FALSE
      )
      
      results <- bind_rows(results, new_row)
      previous_flags <- append(previous_flags, list(new_row))
      
      # Save after each flag
      write_csv(results, output_file)
      
      cat(sprintf("✓ Flagged as %s and saved\n", flag))
      
      i <- i + 1
      
    } else {
      cat("Invalid input. Please try again.\n")
    }
  }
  
  # Final save and summary
  write_csv(results, output_file)
  
  cat("\n=== FLAGGING SUMMARY ===\n")
  cat(sprintf("Output file: %s\n\n", output_file))
  
  flag_summary <- results %>%
    count(flag) %>%
    arrange(flag)
  
  print(flag_summary)
  
  cat("\n")
  
  return(results)
}

#Downloading baseline corrected chlorophyll fluorescence profiles
f <- read_csv(here("files", "8_binAvg-1762199720021.csv"))

#Wrangling CTD profiles, setting date column and renaming columns 
prof <- f %>%
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
         flu_flag = `Fluorometry Chlorophyll flag`)

# Run the flagging tool
# Uncomment the line below when ready to start flagging
qc_results <- flag_fluorescence_profiles(prof)