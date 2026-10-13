# Spring bloom timing and magnitude at KC10
# Input:  files/bloom_timing_kc10.csv (one row per year; add new years there)
# Output: figures/bloom_timing_kc10_no-title.png

library(tidyverse)  # includes readr and lubridate
library(here)

# ---- Data ----
# date  = bloom date, tchla = chlorophyll (mg m-3),
# type  = data source used to identify the bloom that year.
# Any extra columns (e.g. IDs, comments) are read in but not used in the plot.
# Dates are accepted as either YYYY-MM-DD or M/D/YYYY (Excel's default).
sb_timing <- read_csv(here("files", "bloom_timing_kc10_2026-10-07.csv"),
                      col_types = cols(date = col_character(),
                                       id   = col_character())) %>%
  mutate(date = as_date(parse_date_time(date, orders = c("ymd", "mdy"))),
         year = year(date),
         yday = yday(date))

# ---- Plot settings ----
mean_doy   <- round(mean(sb_timing$yday))  # long-term mean bloom day (dashed line)
buoy_start <- 2022                         # first year with daily buoy data
buoy_end   <- 2025                         # last year with daily buoy data

# Years are plotted as discrete levels, so find their y positions
years     <- sort(unique(sb_timing$year))
buoy_pos  <- which(years == buoy_start)
buoy_last <- which(years == buoy_end)

# Size and fill share limits, breaks and title so they merge into one legend
chl_limits <- c(4, 25)
chl_breaks <- scales::pretty_breaks(5)
chl_label  <- expression(atop("Chla", paste("(mg m"^-3, ")")))

# ---- Plot ----
p <- ggplot(sb_timing, aes(x = yday, y = factor(year), size = tchla, fill = tchla)) +
  # Shade the years with daily buoy data
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = buoy_pos - 0.5, ymax = buoy_last + 0.5,
           fill = "grey90", alpha = 0.5) +
  geom_point(pch = 21, colour = "black", stroke = 2) +
  geom_vline(xintercept = mean_doy, linewidth = 1, linetype = "dashed") +
  annotate("text", x = 143, y = (buoy_pos + buoy_last) / 2,  # centred in shading
           label = "Daily Buoy Data", angle = 270, size = 10) +
  scale_size_continuous(range = c(1.5, 28), limits = chl_limits, breaks = chl_breaks) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1,
                       limits = chl_limits, breaks = chl_breaks) +
  guides(fill = guide_legend(), size = guide_legend()) +
  coord_cartesian(xlim = c(40, 150)) +
  labs(x = "YDAY", fill = chl_label, size = chl_label) +
  theme_bw() +
  theme(text = element_text(size = 40),
        axis.title.y = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "right")

ggsave(here("figures", "bloom_timing_kc10_2026-10-07.png"), plot = p,
       width = 9, height = 14, dpi = 300)
