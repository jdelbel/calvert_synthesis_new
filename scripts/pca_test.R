#==============================================================================
# ENVIRONMENTAL PCA ANALYSIS
#==============================================================================
library(tidyverse)
library(here)
library(vegan)
library(patchwork)
library(lubridate)
library(ggrepel)

#==============================================================================
# 1. LOAD DATA
#==============================================================================
env <- read_csv(here("outputs", "enviro_2025-12-19.csv"))
vert <- read_csv(here("outputs", "vertical_structure_metrics.csv"))

#==============================================================================
# 2. DATA PREPARATION
#==============================================================================
env_clean <- env %>% select(-contains("..."))

vert <- vert %>% select(date, int_chla, vertical_concentration_index,
                        scm_intensity)

env_averaged <- env_clean %>%
  group_by(date) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  left_join(vert, by = "date")

pca_data <- env_averaged %>%
  drop_na()

# Add season based on month
pca_data <- pca_data %>%
  mutate(
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall"
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Fall")),
    year = year(date)
  )

cat("\n=== DATA COVERAGE ===\n")
cat("Environmental samples:", nrow(env_clean), "\n")
cat("After joining and dropping NA:", nrow(pca_data), "\n\n")

#==============================================================================
# 3. SELECT ENVIRONMENTAL VARIABLES FOR PCA
#==============================================================================
env_vars <- pca_data %>%
  mutate(sal_diff = sal_30m/sal_1m,
         temp_diff = temp_30m/temp_1m,
         n_diff = nitrate_30m/nitrate_1m) %>% 
  select(
    # Surface conditions
    sal_diff,
    temp_diff,
    # temp_5m,
    # sal_5m,
    # nitrate_1m,
    # nitrate_30m,
    n_diff,
    
    # Stratification
    stratification_index,
    
    
    # Vertical structure
    int_chla,
    vertical_concentration_index,
    scm_intensity,
    
    # Wind/mixing/discharge
    qsm_10,
    wind_b1_3
  ) %>%
  select(where(is.numeric)) %>%
  as.data.frame()

# Check for any remaining NAs or infinite values
cat("\n=== CHECKING FOR MISSING/INFINITE VALUES ===\n")
print(colSums(is.na(env_vars)))
print(colSums(!is.finite(as.matrix(env_vars))))

# Remove any rows with NA or infinite values
env_vars_clean <- env_vars[complete.cases(env_vars) & 
                             apply(env_vars, 1, function(x) all(is.finite(x))), ]

pca_data_clean <- pca_data[complete.cases(env_vars) & 
                             apply(env_vars, 1, function(x) all(is.finite(x))), ]

cat("\n=== FINAL SAMPLE SIZE ===\n")
cat("Samples for PCA:", nrow(env_vars_clean), "\n\n")

#==============================================================================
# 4. RUN PCA
#==============================================================================
# Standardize variables (scale = TRUE centers and scales)
pca_result <- rda(env_vars_clean, scale = TRUE)

cat("\n=== PCA SUMMARY ===\n")
print(summary(pca_result))

# Extract eigenvalues and variance explained
eigenvalues <- pca_result$CA$eig
variance_explained <- eigenvalues / sum(eigenvalues) * 100
cumulative_variance <- cumsum(variance_explained)

cat("\n=== VARIANCE EXPLAINED ===\n")
variance_table <- data.frame(
  PC = paste0("PC", 1:length(eigenvalues)),
  Eigenvalue = eigenvalues,
  Variance_Percent = variance_explained,
  Cumulative_Percent = cumulative_variance
)
print(variance_table)

#==============================================================================
# 5. EXTRACT SCORES FOR PLOTTING
#==============================================================================
# Site scores (samples in PCA space)
site_scores <- scores(pca_result, display = "sites", choices = 1:2)

# Variable loadings (how variables relate to PCs)
var_scores <- scores(pca_result, display = "species", choices = 1:2)

# Prepare plotting data
plot_data <- as.data.frame(site_scores) %>%
  mutate(
    date = pca_data_clean$date,
    season = pca_data_clean$season,
    year = pca_data_clean$year,
    month = pca_data_clean$month
  )

var_arrows <- as.data.frame(var_scores) %>%
  mutate(
    variable = rownames(.),
    label = recode(variable,
                   "temp_5m" = "Temperature",
                   "sal_5m" = "Salinity",
                   "nitrate_5m" = "Nitrate (5m)",
                   "stratification_index" = "Stratification",
                   "int_chla" = "Integrated Chl-a",
                   "vertical_concentration_index" = "Vertical Concentration",
                   "scm_intensity" = "SCM Intensity",
                   "qsm_10" = "Discharge",
                   "wind_b1_3" = "Wind Speed",
                   .default = variable)
  ) %>%
  mutate(
    label_x = PC1 * 1.15,
    label_y = PC2 * 1.15
  )

#==============================================================================
# 6. BIPLOT: PCA WITH ENVIRONMENTAL VECTORS
#==============================================================================
season_colors <- c("Winter" = "#4477AA", "Spring" = "#228833", 
                   "Summer" = "#CCBB44", "Fall" = "#EE6677")

p_pca_season <- ggplot() +
  # Sample points colored by season
  geom_point(data = plot_data,
             aes(x = PC1, y = PC2, fill = season),
             shape = 21, size = 4, alpha = 0.7, color = "black", stroke = 0.5) +
  
  # Variable arrows
  geom_segment(data = var_arrows,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "black", linewidth = 1) +
  
  # Variable labels
  geom_text_repel(data = var_arrows,
                  aes(x = label_x, y = label_y, label = label),
                  color = "black", size = 4, fontface = "bold",
                  box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0) +
  
  scale_fill_manual(values = season_colors, name = "Season") +
  
  labs(
    title = "PCA of Environmental Conditions",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

print(p_pca_season)
ggsave(here("figures", "pca_environmental_season.png"), 
       width = 10, height = 8, dpi = 300, bg = "white")

#==============================================================================
# 7. ALTERNATIVE: COLOR BY YEAR
#==============================================================================
p_pca_year <- ggplot() +
  geom_point(data = plot_data,
             aes(x = PC1, y = PC2, color = as.factor(year)),
             size = 4, alpha = 0.7) +
  
  geom_segment(data = var_arrows,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "black", linewidth = 1) +
  
  geom_text_repel(data = var_arrows,
                  aes(x = label_x, y = label_y, label = label),
                  color = "black", size = 4, fontface = "bold",
                  box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0) +
  
  scale_color_viridis_d(name = "Year") +
  
  labs(
    title = "PCA of Environmental Conditions",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

print(p_pca_year)
ggsave(here("figures", "pca_environmental_year.png"), 
       width = 10, height = 8, dpi = 300, bg = "white")

#==============================================================================
# 8. SCREE PLOT
#==============================================================================
scree_data <- data.frame(
  PC = 1:length(eigenvalues),
  Variance = variance_explained
)

p_scree <- ggplot(scree_data[1:min(10, nrow(scree_data)), ], 
                  aes(x = PC, y = Variance)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Scree Plot",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank()
  )

print(p_scree)
ggsave(here("figures", "pca_scree.png"), 
       width = 8, height = 6, dpi = 300, bg = "white")

#==============================================================================
# 9. HIERARCHICAL CLUSTERING
#==============================================================================

# Use the same standardized environmental data
env_scaled <- scale(env_vars_clean)

# Calculate distance matrix (Euclidean distance on scaled data)
dist_matrix <- dist(env_scaled, method = "euclidean")

# Perform hierarchical clustering (try different methods)
hclust_ward <- hclust(dist_matrix, method = "ward.D2")  # Ward's method
hclust_complete <- hclust(dist_matrix, method = "complete")
hclust_average <- hclust(dist_matrix, method = "average")

# Dendrogram with Ward's method
png(here("figures", "dendrogram_ward.png"), width = 12, height = 8, 
    units = "in", res = 300)
plot(hclust_ward, labels = FALSE, main = "Hierarchical Clustering (Ward's Method)",
     xlab = "Samples", ylab = "Height")
dev.off()

# Determine optimal number of clusters using multiple methods
library(factoextra)

# Elbow method
p_elbow <- fviz_nbclust(env_scaled, FUN = hcut, method = "wss", 
                        hc_method = "ward.D2", k.max = 10) +
  labs(title = "Elbow Method for Optimal Clusters") +
  theme_minimal(base_size = 14)

# Silhouette method
p_silhouette <- fviz_nbclust(env_scaled, FUN = hcut, method = "silhouette",
                             hc_method = "ward.D2", k.max = 10) +
  labs(title = "Silhouette Method for Optimal Clusters") +
  theme_minimal(base_size = 14)

# Gap statistic (can be slow)
p_gap <- fviz_nbclust(env_scaled, FUN = hcut, method = "gap_stat",
                      hc_method = "ward.D2", k.max = 10) +
  labs(title = "Gap Statistic for Optimal Clusters") +
  theme_minimal(base_size = 14)

print(p_gap)

# Combine plots
p_cluster_diagnostics <- p_elbow + p_silhouette
ggsave(here("figures", "cluster_diagnostics.png"), p_cluster_diagnostics,
       width = 12, height = 5, dpi = 300, bg = "white")

# Cut tree into k clusters (adjust k based on diagnostics above)
k <- 4  # Change based on optimal k
clusters <- cutree(hclust_ward, k = k)

# Add clusters to data
pca_data_clean <- pca_data_clean %>%
  mutate(cluster = factor(clusters))

plot_data <- plot_data %>%
  mutate(cluster = factor(clusters))

# PCA biplot colored by cluster
p_pca_cluster <- ggplot() +
  geom_point(data = plot_data,
             aes(x = PC1, y = PC2, fill = cluster),
             shape = 21, size = 4, alpha = 0.7, color = "black", stroke = 0.5) +
  
  geom_segment(data = var_arrows,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "black", linewidth = 1) +
  
  geom_text_repel(data = var_arrows,
                  aes(x = label_x, y = label_y, label = label),
                  color = "black", size = 4, fontface = "bold",
                  box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0) +
  
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  
  labs(
    title = "PCA with Hierarchical Clusters",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

print(p_pca_cluster)
ggsave(here("figures", "pca_clusters.png"), 
       width = 10, height = 8, dpi = 300, bg = "white")

# Characterize clusters by environmental variables
cluster_summary <- pca_data_clean %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    across(c(temp_5m, sal_5m, nitrate_1m, stratification_index, 
             int_chla, vertical_concentration_index, scm_intensity,
             qsm_10, wind_b1_3),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE))),
    .groups = 'drop'
  )

print(cluster_summary)
write_csv(cluster_summary, here("outputs", "cluster_environmental_summary.csv"))

# Seasonal composition of clusters
cluster_season <- pca_data_clean %>%
  count(cluster, season) %>%
  group_by(cluster) %>%
  mutate(proportion = n / sum(n))

p_cluster_season <- ggplot(cluster_season, 
                           aes(x = cluster, y = proportion, fill = season)) +
  geom_col(position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = season_colors, name = "Season") +
  labs(
    title = "Seasonal Composition of Environmental Clusters",
    x = "Cluster",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(p_cluster_season)
ggsave(here("figures", "cluster_seasonal_composition.png"),
       width = 8, height = 6, dpi = 300, bg = "white")

#==============================================================================
# 10. BOXPLOTS OF MAIN PCA DRIVERS BY CLUSTER
#==============================================================================

# Identify top loadings for PC1 and PC2
pc1_loadings <- var_scores[, "PC1"]
pc2_loadings <- var_scores[, "PC2"]

cat("\n=== TOP PC1 LOADINGS ===\n")
print(sort(abs(pc1_loadings), decreasing = TRUE))

cat("\n=== TOP PC2 LOADINGS ===\n")
print(sort(abs(pc2_loadings), decreasing = TRUE))

# Prepare data for plotting - combine cluster assignments with environmental data
cluster_env_data <- pca_data_clean %>%
  select(cluster, temp_5m, sal_5m, nitrate_1m, nitrate_10m, nitrate_30m,
         stratification_index, int_chla, vertical_concentration_index, 
         scm_intensity, qsm_10, wind_b1_3)

# Pivot to long format for faceted plotting
cluster_env_long <- cluster_env_data %>%
  pivot_longer(cols = -cluster, names_to = "variable", values_to = "value") %>%
  mutate(
    variable_label = recode(variable,
                            "temp_5m" = "Temperature (°C)",
                            "sal_5m" = "Salinity (PSU)",
                            "nitrate_1m" = "Nitrate 1m (µM)",
                            "nitrate_10m" = "Nitrate 10m (µM)",
                            "nitrate_30m" = "Nitrate 30m (µM)",
                            "stratification_index" = "Stratification (kg/m³)",
                            "int_chla" = "Integrated Chl-a (mg/m²)",
                            "vertical_concentration_index" = "Vertical Conc. Index",
                            "scm_intensity" = "SCM Intensity",
                            "qsm_10" = "Discharge (m³/s)",
                            "wind_b1_3" = "Wind Speed (m/s)")
  )

# Select key drivers based on PCA loadings (adjust based on your results)
key_drivers <- c("temp_5m", "sal_5m", "nitrate_1m", "stratification_index", 
                 "int_chla", "scm_intensity")

cluster_env_key <- cluster_env_long %>%
  filter(variable %in% key_drivers) %>%
  mutate(variable_label = factor(variable_label, 
                                 levels = c("Temperature (°C)", 
                                            "Salinity (PSU)",
                                            "Nitrate 1m (µM)",
                                            "Stratification (kg/m³)",
                                            "Integrated Chl-a (mg/m²)",
                                            "SCM Intensity")))

# Faceted boxplots
p_boxplots_facet <- ggplot(cluster_env_key, 
                           aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.alpha = 0.5) +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 3) +
  labs(
    title = "Environmental Drivers by Cluster",
    x = "Cluster",
    y = "Value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none",
    axis.text.x = element_text(size = 10)
  )

print(p_boxplots_facet)
ggsave(here("figures", "cluster_boxplots_faceted.png"), 
       width = 12, height = 8, dpi = 300, bg = "white")

# Alternative: Individual plots for top 4 drivers (better for presentations)
# Temperature
p_temp <- ggplot(cluster_env_data, aes(x = cluster, y = temp_5m, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Temperature", y = "Temperature (°C)", x = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

# Nitrate
p_nitrate <- ggplot(cluster_env_data, aes(x = cluster, y = nitrate_1m, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Surface Nitrate", y = "Nitrate 1m (µM)", x = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

# Integrated chlorophyll
p_chla <- ggplot(cluster_env_data, aes(x = cluster, y = int_chla, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Integrated Chlorophyll-a", y = "Integrated Chl-a (mg/m²)", x = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

# Stratification
p_strat <- ggplot(cluster_env_data, aes(x = cluster, y = stratification_index, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Stratification", y = "Stratification Index (kg/m³)", x = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

# Combine into 2x2 panel
p_boxplots_panel <- (p_temp + p_nitrate) / (p_chla + p_strat)

print(p_boxplots_panel)
ggsave(here("figures", "cluster_boxplots_panel.png"), 
       width = 10, height = 8, dpi = 300, bg = "white")

# Statistical comparison between clusters (Kruskal-Wallis tests)
kw_results <- cluster_env_data %>%
  select(-cluster) %>%
  map_df(~{
    test <- kruskal.test(.x ~ cluster_env_data$cluster)
    tibble(
      p_value = test$p.value,
      statistic = test$statistic
    )
  }, .id = "variable") %>%
  mutate(
    significant = p_value < 0.05,
    variable_label = recode(variable,
                            "temp_5m" = "Temperature",
                            "sal_5m" = "Salinity",
                            "nitrate_1m" = "Nitrate 1m",
                            "nitrate_10m" = "Nitrate 10m",
                            "nitrate_30m" = "Nitrate 30m",
                            "stratification_index" = "Stratification",
                            "int_chla" = "Integrated Chl-a",
                            "vertical_concentration_index" = "Vertical Conc. Index",
                            "scm_intensity" = "SCM Intensity",
                            "qsm_10" = "Discharge",
                            "wind_b1_3" = "Wind Speed")
  ) %>%
  arrange(p_value)

cat("\n=== KRUSKAL-WALLIS TESTS FOR CLUSTER DIFFERENCES ===\n")
print(kw_results)

write_csv(kw_results, here("outputs", "cluster_kruskal_wallis_tests.csv"))

# Faceted boxplots - ALL VARIABLES
p_boxplots_facet_all <- ggplot(cluster_env_long, 
                               aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.alpha = 0.5) +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 3) +
  labs(
    title = "Environmental Drivers by Cluster",
    x = "Cluster",
    y = "Value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none",
    axis.text.x = element_text(size = 10)
  )

print(p_boxplots_facet_all)
ggsave(here("figures", "cluster_boxplots_all_variables.png"), 
       width = 14, height = 12, dpi = 300, bg = "white")

# Year composition of clusters
cluster_year <- pca_data_clean %>%
  count(cluster, year) %>%
  group_by(cluster) %>%
  mutate(proportion = n / sum(n))

# Stacked bar plot by cluster
p_cluster_year <- ggplot(cluster_year, 
                         aes(x = cluster, y = proportion, fill = factor(year))) +
  geom_col(position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_viridis_d(name = "Year", option = "turbo") +
  labs(
    title = "Interannual Composition of Environmental Clusters",
    x = "Cluster",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(p_cluster_year)
ggsave(here("figures", "cluster_year_composition.png"),
       width = 8, height = 6, dpi = 300, bg = "white")

# Alternative: Year on x-axis to see temporal trends
cluster_year_trend <- pca_data_clean %>%
  count(year, cluster) %>%
  group_by(year) %>%
  mutate(proportion = n / sum(n))

p_year_cluster_trend <- ggplot(cluster_year_trend, 
                               aes(x = year, y = proportion, fill = cluster)) +
  geom_col(position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  scale_x_continuous(breaks = unique(cluster_year_trend$year)) +
  labs(
    title = "Temporal Trends in Environmental Clusters",
    x = "Year",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_year_cluster_trend)
ggsave(here("figures", "year_cluster_trend.png"),
       width = 10, height = 6, dpi = 300, bg = "white")

# Combined seasonal and interannual view
cluster_year_season <- pca_data_clean %>%
  count(cluster, year, season) %>%
  group_by(cluster, year) %>%
  mutate(proportion = n / sum(n))

p_cluster_year_season <- ggplot(cluster_year_season,
                                aes(x = year, y = proportion, fill = season)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = season_colors, name = "Season") +
  scale_x_continuous(breaks = unique(cluster_year_season$year)) +
  facet_wrap(~ cluster, ncol = 2) +
  labs(
    title = "Seasonal and Interannual Patterns by Cluster",
    x = "Year",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

print(p_cluster_year_season)
ggsave(here("figures", "cluster_year_season_composition.png"),
       width = 12, height = 8, dpi = 300, bg = "white")

# Summary table
cluster_year_summary <- pca_data_clean %>%
  count(cluster, year) %>%
  pivot_wider(names_from = year, values_from = n, values_fill = 0)

print(cluster_year_summary)
write_csv(cluster_year_summary, here("outputs", "cluster_year_summary.csv"))
