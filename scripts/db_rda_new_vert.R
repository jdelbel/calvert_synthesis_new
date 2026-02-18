#==============================================================================
# CLEAN db-RDA WITH PROPER VARIABLE SELECTION
#==============================================================================

library(tidyverse)
library(here)
library(vegan)
library(patchwork)
library(lubridate)

#==============================================================================
# 1. LOAD DATA
#==============================================================================

chem <- read_csv(here("outputs", "chem_clusters_2025-12-03.csv"))
env <- read_csv(here("outputs", "enviro_2025-12-19.csv"))
vert <- read_csv(here("outputs", "vert_fl_metrics_2026-01-14.csv"))

#==============================================================================
# 2. DATA PREPARATION
#==============================================================================

chem_clean <- chem %>% select(-contains("..."))
env_clean <- env %>% select(-contains("..."))

vert <- vert %>% select(date, int_chla = f_integrated,
                        peakiness,
                        max = f_max,
                        depth_max = depth_max)

env_averaged <- env_clean %>%
  group_by(date) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  left_join(vert, by = "date")

rda_data <- chem_clean %>%
  left_join(env_averaged, by = "date") %>%
  drop_na()

# Add season based on month
rda_data <- rda_data %>%
  mutate(
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall"
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Fall"))
  )

cat("\n=== DATA COVERAGE ===\n")
cat("CHEMTAX samples:", nrow(chem_clean), "\n")
cat("After joining:", nrow(rda_data), "\n\n")

#==============================================================================
# 3. SELECT VARIABLES
#==============================================================================

phyto_groups <- c("diat", "raph", "dict", "dino", "cryp", "GA", "hapt", "cyan")
phyto_response <- rda_data %>% select(all_of(phyto_groups)) %>% as.data.frame()

env_predictors <- rda_data %>%
  select(
    # Surface conditions
    temp_5m,
    sal_5m,
    nitrate_5m,
    # nitrate_10m,
    
    # Vertical gradients
    # delta_temp_5_30,
    # delta_sigma_2_30,
    # mean_N2,
    stratification_index,
    
    # Vertical structure
    int_chla,
    # centroid_depth,
    # profile_shape,
    # surface_depletion,
    # scm_depth,
    # vertical_gradient,
    peakiness,
    # max,
    depth_max,
    
    # Wind/mixing
    qsm_10,
    wind_b1_3,
    wd_b1_3
  ) %>%
  select(where(is.numeric)) %>%
  as.data.frame()

# Remove zero-variance columns
zero_var <- env_predictors %>%
  summarise(across(everything(), ~var(.x, na.rm = TRUE) == 0 | is.na(var(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), values_to = "is_zero") %>%
  filter(is_zero) %>%
  pull(name)

if(length(zero_var) > 0) {
  cat("Removing zero-variance:", paste(zero_var, collapse = ", "), "\n")
  env_predictors <- env_predictors %>% select(-all_of(zero_var))
}

env_scaled <- scale(env_predictors)

cat("=== VARIABLES AVAILABLE ===\n")
cat(paste(colnames(env_scaled), collapse = ", "), "\n\n")

#==============================================================================
# 4. RUN db-RDA WITH STRICTER SELECTION
#==============================================================================

phyto_dist <- vegdist(phyto_response, method = "bray")

# Full model
dbrda_full <- capscale(phyto_dist ~ ., data = as.data.frame(env_scaled))

cat("=== VIF - FULL MODEL ===\n")
vif_full <- vif.cca(dbrda_full)
print(sort(vif_full, decreasing = TRUE))
cat("\n")

# Forward selection with STRICTER p-value threshold
dbrda_null <- capscale(phyto_dist ~ 1, data = as.data.frame(env_scaled))

dbrda_forward <- ordiR2step(
  dbrda_null,
  scope = formula(dbrda_full),
  direction = "forward",
  permutations = 9999,
  Pin = 0.05,  # Only add variables with p < 0.05
  trace = TRUE
)

cat("\n=== SELECTED VARIABLES ===\n")
selected_vars <- attr(terms(dbrda_forward), "term.labels")
cat(paste(selected_vars, collapse = ", "), "\n\n")

# Create final model with ONLY significant variables (p < 0.05)
cat("=== TESTING EACH VARIABLE INDIVIDUALLY ===\n")
anova_terms <- anova(dbrda_forward, by = "terms", permutations = 999)
print(anova_terms)

# Keep only variables with p < 0.05
sig_vars <- rownames(anova_terms)[anova_terms$`Pr(>F)` < 0.05]
sig_vars <- sig_vars[!is.na(sig_vars) & sig_vars != "Residual"]

cat("\n=== SIGNIFICANT VARIABLES (p < 0.05) ===\n")
cat(paste(sig_vars, collapse = ", "), "\n\n")

# Rebuild final model with only significant variables
if(length(sig_vars) > 0) {
  formula_final <- as.formula(paste("phyto_dist ~", paste(sig_vars, collapse = " + ")))
  dbrda_final <- capscale(formula_final, data = as.data.frame(env_scaled))
} else {
  cat("WARNING: No significant variables!\n")
  dbrda_final <- dbrda_forward
}

#==============================================================================
# 5. FINAL MODEL STATISTICS
#==============================================================================

r2 <- RsquareAdj(dbrda_final)
variance <- summary(dbrda_final)$cont$importance

cat("=== FINAL MODEL ===\n")
cat("Variables:", paste(attr(terms(dbrda_final), "term.labels"), collapse = ", "), "\n")
cat("Adjusted R² =", round(r2$adj.r.squared, 3), "\n")
cat("CAP1:", round(variance[2, 1] * 100, 1), "%\n")
cat("CAP2:", round(variance[2, 2] * 100, 1), "%\n\n")

# Model significance
cat("=== MODEL SIGNIFICANCE ===\n")
print(anova(dbrda_final, permutations = 9999))

# Axis significance
cat("\n=== AXIS SIGNIFICANCE ===\n")
print(anova(dbrda_final, by = "axis", permutations = 9999))

# Variable significance
cat("\n=== VARIABLE SIGNIFICANCE ===\n")
print(anova(dbrda_final, by = "terms", permutations = 9999))

# VIF for final model
cat("\n=== VIF - FINAL MODEL ===\n")
vif_final <- vif.cca(dbrda_final)
vif_df <- data.frame(
  Variable = names(vif_final),
  VIF = vif_final
) %>% arrange(desc(VIF))
print(vif_df)

#==============================================================================
# 5B. VARIABLE IMPORTANCE DIAGNOSTICS
#==============================================================================

cat("\n=== DETAILED VARIABLE IMPORTANCE DIAGNOSTICS ===\n\n")

# 1. Biplot scores (arrow lengths = importance)
cat("1. BIPLOT SCORES (arrow length indicates importance):\n")
biplot_scores <- scores(dbrda_final, display = "bp", choices = 1:2)
biplot_lengths <- sqrt(biplot_scores[,1]^2 + biplot_scores[,2]^2)
biplot_df <- data.frame(
  Variable = rownames(biplot_scores),
  CAP1 = biplot_scores[,1],
  CAP2 = biplot_scores[,2],
  Length = biplot_lengths
) %>% arrange(desc(Length))
print(biplot_df)

# 2. Correlations with ordination axes
cat("\n2. CORRELATIONS WITH ORDINATION AXES:\n")
site_scores_df <- scores(dbrda_final, display = "sites", choices = 1:2)
cors_with_axes <- cor(env_scaled[, sig_vars], site_scores_df)
cat("\nCAP1 correlations:\n")
print(sort(cors_with_axes[,1]^2, decreasing = TRUE))  # R-squared
cat("\nCAP2 correlations:\n")
print(sort(cors_with_axes[,2]^2, decreasing = TRUE))  # R-squared

# 3. Variance partitioning - each variable's unique contribution
cat("\n3. VARIANCE PARTITIONING (unique contribution of each variable):\n")
var_importance <- data.frame(
  Variable = character(),
  R2_full = numeric(),
  R2_without = numeric(),
  Unique_R2 = numeric(),
  Percent_of_total = numeric()
)

r2_full <- RsquareAdj(dbrda_final)$r.squared

for(var in sig_vars) {
  # Fit model without this variable
  other_vars <- sig_vars[sig_vars != var]
  if(length(other_vars) > 0) {
    formula_reduced <- as.formula(paste("phyto_dist ~", paste(other_vars, collapse = " + ")))
    dbrda_reduced <- capscale(formula_reduced, data = as.data.frame(env_scaled))
    r2_without <- RsquareAdj(dbrda_reduced)$r.squared
  } else {
    r2_without <- 0
  }
  
  unique_r2 <- r2_full - r2_without
  var_importance <- rbind(var_importance, data.frame(
    Variable = var,
    R2_full = r2_full,
    R2_without = r2_without,
    Unique_R2 = unique_r2,
    Percent_of_total = (unique_r2 / r2_full) * 100
  ))
}

var_importance <- var_importance %>% arrange(desc(Unique_R2))
print(var_importance)

# 4. Marginal and conditional effects
cat("\n4. MARGINAL vs CONDITIONAL EFFECTS:\n")
cat("(Marginal = variable alone; Conditional = after other variables)\n\n")

marginal_r2 <- numeric()
for(var in sig_vars) {
  formula_single <- as.formula(paste("phyto_dist ~", var))
  dbrda_single <- capscale(formula_single, data = as.data.frame(env_scaled))
  marginal_r2[var] <- RsquareAdj(dbrda_single)$r.squared
}

effects_df <- data.frame(
  Variable = sig_vars,
  Marginal_R2 = marginal_r2[sig_vars],
  Unique_R2 = var_importance$Unique_R2,
  Shared_R2 = marginal_r2[sig_vars] - var_importance$Unique_R2
) %>% arrange(desc(Marginal_R2))

print(effects_df)

cat("\n=== INTERPRETATION ===")
cat("\n- Marginal R²: How much variance the variable explains alone")
cat("\n- Unique R²: How much it adds beyond other variables")
cat("\n- Shared R²: Variance explained jointly with other variables")
cat("\n- High Shared/Low Unique = redundant with other variables")
cat("\n- High Unique = important independent driver\n\n")

# Cluster separation
cat("\n=== CLUSTER SEPARATION ===\n")
site_scores <- scores(dbrda_final, display = "sites", choices = 1:2)
cluster_cap <- data.frame(
  CAP1 = site_scores[, 1],
  CAP2 = site_scores[, 2],
  cluster = factor(rda_data$abs_clust)
)

cat("\nCAP1 ~ Cluster:\n")
print(summary(aov(CAP1 ~ cluster, data = cluster_cap)))

cat("\nCAP2 ~ Cluster:\n")
print(summary(aov(CAP2 ~ cluster, data = cluster_cap)))

#==============================================================================
# 6. PREPARE SPECIES VECTORS FOR PLOTTING
#==============================================================================

# For db-RDA, use envfit to properly fit species vectors
species_fit <- envfit(dbrda_final, phyto_response, permutations = 999)

cat("\n=== SPECIES VECTOR FITTING ===\n")
print(species_fit)

# Extract the scores
species_arrows <- as.data.frame(scores(species_fit, display = "vectors")) %>%
  rownames_to_column("group") %>%
  rename(CAP1 = CAP1, CAP2 = CAP2) %>%
  mutate(
    group_label = recode(group,
                         "diat" = "Diatoms",
                         "raph" = "Raphidophytes",
                         "dict" = "Dictyochophytes",
                         "dino" = "Dinoflagellates",
                         "cryp" = "Cryptophytes",
                         "GA" = "Green Algae",
                         "hapt" = "Haptophytes",
                         "cyan" = "Cyanobacteria")
  )

# Get p-values and R² for reference
species_stats <- data.frame(
  group = names(species_fit$vectors$r),
  r2 = species_fit$vectors$r,
  pval = species_fit$vectors$pvals
)

species_arrows <- species_arrows %>%
  left_join(species_stats, by = "group")

cat("\n=== SPECIES VECTORS (all groups) ===\n")
print(species_arrows)

#==============================================================================
# 6B. FILTERED SPECIES VECTORS WITH GGREPEL
#==============================================================================

library(ggrepel)

# Filter to significant species only (adjust threshold as needed)
species_arrows_filtered <- species_arrows %>%
  filter(pval < 0.05)  # You can also add: | r2 > 0.3

cat("\n=== FILTERED SPECIES VECTORS (p < 0.05) ===\n")
print(species_arrows_filtered)

cluster_colors <- c("1" = "#DDAA33", "2" = "#228833", "3" = "#CC6677", "4" = "#4477AA")
cluster_labels <- c("1" = "1 (Mixed)", "2" = "2 (Diatoms)", 
                    "3" = "3 (Dinoflagellates)", "4" = "4 (Winter)")

plot_data <- as.data.frame(site_scores) %>%
  mutate(
    cluster = factor(rda_data$abs_clust),
    season = rda_data$season,
    scm_intensity = rda_data$scm_intensity
  )

env_arrows <- as.data.frame(scores(dbrda_final, display = "bp", choices = 1:2)) %>%
  mutate(
    variable = rownames(.),
    label = recode(variable,
                   "nitrate_5m" = "Nitrate (5m)",
                   "temp_5m" = "Temperature",
                   "sal_5m" = "Salinity",
                   "int_chla" = "Integrated Chl-a",
                   "scm_depth" = "SCM Depth",
                   "scm_intensity" = "SCM Intensity",
                   "delta_temp_5_30" = "ΔTemp",
                   "delta_sigma_2_30" = "ΔDensity",
                   "qsm_10" = "Snow-Mountain Discharge",
                   "vertical_gradient" = "Vertical Gradient",
                   .default = variable)
  ) %>%
  mutate(
    label_x = CAP1 * 4.2,
    label_y = CAP2 * 4.2
  )

var_exp <- summary(dbrda_final)$cont$importance[2, 1:2] * 100

p_dbrda <- ggplot() +
  geom_point(data = plot_data,
             aes(x = CAP1, y = CAP2, fill = cluster, shape = cluster),
             size = 7, alpha = 0.85, color = "black", stroke = 1) +
  
  # Environmental arrows (black)
  geom_segment(data = env_arrows,
               aes(x = 0, y = 0, xend = CAP1 * 4, yend = CAP2 * 4),
               arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
               color = "red", linewidth = 1.8) +
  
  # Environmental labels with ggrepel
  geom_text_repel(data = env_arrows,
                  aes(x = label_x, y = label_y, label = label),
                  color = "red", size = 6, fontface = "bold",
                  box.padding = 0.5, 
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  segment.size = 0) +
  
  # Species arrows (blue, dashed) - FILTERED and scaled at 3.5
  geom_segment(data = species_arrows_filtered %>% 
                 mutate(arrow_x = CAP1 * 3.5, arrow_y = CAP2 * 3.5),
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "#0072B2", linewidth = 1.8, linetype = "dashed", alpha = 0.7) +
  
  # Species labels with ggrepel
  geom_text_repel(data = species_arrows_filtered %>% 
                    mutate(label_x = CAP1 * 3.7, label_y = CAP2 * 3.7),
                  aes(x = label_x, y = label_y, label = group_label),
                  color = "#0072B2", size = 8, fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "#0072B2",
                  segment.size = 0,
                  segment.linetype = "dotted") +
  
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
                     labels = cluster_labels, name = "Community Cluster") +
  scale_fill_manual(values = cluster_colors, labels = cluster_labels,
                    name = "Community Cluster") +
  
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  
  labs(
    title = "Environmental Drivers of Phytoplankton Community Composition",
    x = paste0("CAP1 (", round(var_exp[1], 1), "%)"),
    y = paste0("CAP2 (", round(var_exp[2], 1), "%)"),
    caption = paste0("db-RDA (Bray-Curtis) | R² = ", round(r2$adj.r.squared, 3),
                     " | ", length(sig_vars), " significant variables | Showing ", 
                     nrow(species_arrows_filtered), " significant species vectors (p < 0.05)")
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(25, 25, 25, 25),
    plot.clip = "off"
  ) +
  coord_cartesian(clip = "off")

print(p_dbrda)
ggsave(here("figures", "dbrda_final_new_vert.png"),
       width = 14, height = 10, dpi = 400, bg = "white")

#==============================================================================
# 6C. VARIABLE IMPORTANCE VISUALIZATIONS
#==============================================================================

# Create importance plot
p_importance <- ggplot(var_importance, aes(x = reorder(Variable, Unique_R2), y = Unique_R2 * 100)) +
  geom_col(fill = "#4477AA", alpha = 0.8) +
  geom_text(aes(label = paste0(round(Unique_R2 * 100, 1), "%")), 
            hjust = -0.2, size = 5, fontface = "bold") +
  coord_flip() +
  labs(
    title = "Unique Contribution of Each Variable",
    subtitle = "Variance explained beyond other variables",
    x = NULL,
    y = "Unique R² (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Marginal vs Unique comparison
p_marginal_unique <- effects_df %>%
  pivot_longer(cols = c(Marginal_R2, Unique_R2), names_to = "Type", values_to = "R2") %>%
  mutate(Type = recode(Type, "Marginal_R2" = "Marginal", "Unique_R2" = "Unique")) %>%
  ggplot(aes(x = reorder(Variable, R2), y = R2 * 100, fill = Type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Marginal" = "#DDAA33", "Unique" = "#4477AA")) +
  coord_flip() +
  labs(
    title = "Marginal vs Unique Effects",
    subtitle = "High marginal but low unique = redundant with other variables",
    x = NULL,
    y = "R² (%)",
    fill = "Effect Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save importance figures
ggsave(here("figures", "variable_importance_unique_new_vert.png"), p_importance, 
       width = 10, height = 6, dpi = 300, bg = "white")
ggsave(here("figures", "variable_importance_comparison_new_vert.png"), p_marginal_unique, 
       width = 10, height = 6, dpi = 300, bg = "white")

#==============================================================================
# 7. THREE-PANEL FIGURE WITH SPECIES VECTORS
#==============================================================================

# Define color schemes
season_colors <- c(
  "Winter" = "#4575B4",   # Blue
  "Spring" = "#7FBC41",   # Green
  "Summer" = "#F46D43",   # Orange-red
  "Fall" = "#A6611A"      # Brown
)

# Calculate appropriate arrow scaling to keep vectors inside plot
arrow_scale <- 1.2

# Recalculate label positions with new scaling
env_arrows_panel <- env_arrows %>%
  mutate(
    arrow_x = CAP1 * arrow_scale,
    arrow_y = CAP2 * arrow_scale,
    label_x = CAP1 * (arrow_scale + 0.15),
    label_y = CAP2 * (arrow_scale + 0.15)
  )

# Common theme for all panels
panel_theme <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.box.spacing = unit(0.3, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    plot.margin = margin(10, 15, 10, 15)
  )

# Panel A: Clusters
p_clusters <- ggplot() +
  geom_point(data = plot_data,
             aes(x = CAP1, y = CAP2, fill = cluster, shape = cluster),
             size = 5.5, alpha = 0.85, color = "black", stroke = 0.8) +
  
  # Environmental arrows (black)
  geom_segment(data = env_arrows_panel,
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
               color = "black", linewidth = 1.2) +
  
  geom_text(data = env_arrows_panel,
            aes(x = label_x, y = label_y, label = label),
            color = "black", size = 4.5, fontface = "bold") +
  
  # Species arrows (blue, dashed) - scaled to fit within bounds
  geom_segment(data = species_arrows %>% mutate(arrow_x = CAP1 * 0.4, arrow_y = CAP2 * 0.4),
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "#0072B2", linewidth = 1, linetype = "dashed", alpha = 0.7) +
  
  geom_text(data = species_arrows %>% mutate(label_x = CAP1 * 0.5, label_y = CAP2 * 0.5),
            aes(x = label_x, y = label_y, label = group_label),
            color = "#0072B2", size = 3.5, fontface = "italic") +
  
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
                     labels = cluster_labels, name = "Community Cluster") +
  scale_fill_manual(values = cluster_colors, labels = cluster_labels,
                    name = "Community Cluster") +
  
  labs(
    title = "(A) Community Clusters",
    x = paste0("CAP1 (", round(var_exp[1], 1), "%)"),
    y = paste0("CAP2 (", round(var_exp[2], 1), "%)")
  ) +
  
  guides(
    fill = guide_legend(nrow = 1, override.aes = list(size = 5)),
    shape = guide_legend(nrow = 1, override.aes = list(size = 5))
  ) +
  
  panel_theme

# Panel B: Season (with cluster shapes retained)
p_season <- ggplot() +
  geom_point(data = plot_data,
             aes(x = CAP1, y = CAP2, fill = season, shape = cluster),
             size = 5.5, alpha = 0.85, color = "black", stroke = 0.8) +
  
  # Environmental arrows (black)
  geom_segment(data = env_arrows_panel,
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
               color = "black", linewidth = 1.2) +
  
  geom_text(data = env_arrows_panel,
            aes(x = label_x, y = label_y, label = label),
            color = "black", size = 4.5, fontface = "bold") +
  
  # Species arrows (blue, dashed) - scaled to fit within bounds
  geom_segment(data = species_arrows %>% mutate(arrow_x = CAP1 * 0.4, arrow_y = CAP2 * 0.4),
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "#0072B2", linewidth = 1, linetype = "dashed", alpha = 0.7) +
  
  geom_text(data = species_arrows %>% mutate(label_x = CAP1 * 0.5, label_y = CAP2 * 0.5),
            aes(x = label_x, y = label_y, label = group_label),
            color = "#0072B2", size = 3.5, fontface = "italic") +
  
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
                     guide = "none") +  # Hide shape legend
  scale_fill_manual(values = season_colors, name = "Season") +
  
  labs(
    title = "(B) Seasonal Pattern",
    x = paste0("CAP1 (", round(var_exp[1], 1), "%)"),
    y = ""  # Remove y-axis label
  ) +
  
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 5, shape = 21))) +
  
  panel_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Panel C: SCM Intensity (with cluster shapes retained)
p_scm <- ggplot() +
  geom_point(data = plot_data,
             aes(x = CAP1, y = CAP2, fill = scm_intensity, shape = cluster),
             size = 5.5, alpha = 0.85, color = "black", stroke = 0.8) +
  
  # Environmental arrows (black)
  geom_segment(data = env_arrows_panel,
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
               color = "black", linewidth = 1.2) +
  
  geom_text(data = env_arrows_panel,
            aes(x = label_x, y = label_y, label = label),
            color = "black", size = 4.5, fontface = "bold") +
  
  # Species arrows (blue, dashed) - scaled to fit within bounds
  geom_segment(data = species_arrows %>% mutate(arrow_x = CAP1 * 0.4, arrow_y = CAP2 * 0.4),
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "#0072B2", linewidth = 1, linetype = "dashed", alpha = 0.7) +
  
  geom_text(data = species_arrows %>% mutate(label_x = CAP1 * 0.5, label_y = CAP2 * 0.5),
            aes(x = label_x, y = label_y, label = group_label),
            color = "#0072B2", size = 3.5, fontface = "italic") +
  
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
                     guide = "none") +  # Hide shape legend
  scale_fill_viridis_c(name = "SCM Intensity", option = "plasma") +
  
  labs(
    title = "(C) Subsurface Chlorophyll Maximum",
    x = paste0("CAP1 (", round(var_exp[1], 1), "%)"),
    y = ""  # Remove y-axis label
  ) +
  
  guides(fill = guide_colorbar(barwidth = 15, barheight = 1, title.position = "top",
                               title.hjust = 0.5)) +
  
  panel_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Combine panels
p_combined <- p_clusters + p_season + p_scm + 
  plot_layout(ncol = 3, widths = c(1, 1, 1)) +
  plot_annotation(
    title = "Environmental Drivers of Phytoplankton Community Composition",
    subtitle = paste0("db-RDA (Bray-Curtis) | Adjusted R² = ", round(r2$adj.r.squared, 3),
                      " | ", length(sig_vars), " significant variables"),
    theme = theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 15))
    )
  )

print(p_combined)
ggsave(here("figures", "dbrda_3panel_new_vert.png"), width = 24, height = 10, dpi = 400, bg = "white")

#==============================================================================
# 8. INVESTIGATE CLUSTER 2 HETEROGENEITY
#==============================================================================

# Add ordination coordinates to your data
plot_data_full <- plot_data %>%
  mutate(
    date = rda_data$date,
    diat_prop = rda_data$diat / rowSums(rda_data[, phyto_groups]),
    temp_5m = rda_data$temp_5m,
    int_chla = rda_data$int_chla,
    scm_intensity = rda_data$scm_intensity,
    nitrate_5m = rda_data$nitrate_5m
  )

# Identify NE vs SE diatom samples (Cluster 2 only)
cluster2_samples <- plot_data_full %>%
  filter(cluster == "2") %>%
  mutate(
    quadrant = case_when(
      CAP1 > 0 & CAP2 > 0 ~ "NE (High Chl-a)",
      CAP1 > 0 & CAP2 < 0 ~ "SE (Warm/SCM)",
      TRUE ~ "Other"
    )
  ) %>%
  filter(quadrant != "Other")

cat("\n=== CLUSTER 2 SPLIT BY QUADRANT ===\n")
table(cluster2_samples$quadrant)

# Compare environmental conditions
cat("\n=== ENVIRONMENTAL CONDITIONS - CLUSTER 2 ===\n")
cluster2_summary <- cluster2_samples %>%
  group_by(quadrant) %>%
  summarise(
    n = n(),
    mean_temp = mean(temp_5m, na.rm = TRUE),
    mean_chla = mean(int_chla, na.rm = TRUE),
    mean_scm = mean(scm_intensity, na.rm = TRUE),
    mean_nitrate = mean(nitrate_5m, na.rm = TRUE),
    mean_diat_prop = mean(diat_prop, na.rm = TRUE),
    .groups = 'drop'
  )
print(cluster2_summary)

# Check temporal patterns
cat("\n=== TEMPORAL DISTRIBUTION ===\n")
cluster2_temporal <- cluster2_samples %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(quadrant) %>%
  summarise(
    date_range = paste(min(date), "to", max(date)),
    mean_month = mean(month),
    .groups = 'drop'
  )
print(cluster2_temporal)

# Visualize the split
p_cluster2_split <- ggplot(cluster2_samples, aes(x = int_chla, y = temp_5m, color = quadrant)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = c("NE (High Chl-a)" = "#228833", 
                                "SE (Warm/SCM)" = "red")) +
  labs(
    title = "Cluster 2 (Diatoms): Environmental Conditions",
    x = "Integrated Chlorophyll-a (mg/m²)",
    y = "Temperature (5m, °C)",
    color = "Ordination Position"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

print(p_cluster2_split)
ggsave(here("figures", "cluster2_split_investigation.png"), 
       width = 10, height = 7, dpi = 300, bg = "white")

# Check community composition differences
cat("\n=== PHYTOPLANKTON COMPOSITION BY QUADRANT ===\n")
composition_comparison <- cluster2_samples %>%
  left_join(rda_data %>% select(date, all_of(phyto_groups)), by = "date") %>%
  group_by(quadrant) %>%
  summarise(across(all_of(phyto_groups), ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  mutate(across(all_of(phyto_groups), ~.x / rowSums(across(all_of(phyto_groups))) * 100))

print(composition_comparison)

# Time series of ordination position
p_cluster2_time <- ggplot(cluster2_samples, aes(x = date, y = CAP2, color = quadrant)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("NE (High Chl-a)" = "#228833", 
                                "SE (Warm/SCM)" = "#117733")) +
  labs(
    title = "Cluster 2: Temporal Pattern of Ordination Position",
    x = "Date",
    y = "CAP2 Position",
    color = "Quadrant"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

print(p_cluster2_time)
ggsave(here("figures", "cluster2_temporal_pattern_new_vert.png"), 
       width = 12, height = 6, dpi = 300, bg = "white")

cat("\n✅ Analysis complete!\n")
cat("Final model has", length(sig_vars), "significant variables (p < 0.05)\n")
cat("\n=== OUTPUT FILES ===\n")
cat("Individual plot: dbrda_final.png\n")
cat("3-panel figure with species vectors: dbrda_3panel.png\n")
cat("Variable importance (unique contributions): variable_importance_unique.png\n")
cat("Variable importance (marginal vs unique): variable_importance_comparison.png\n")
cat("Cluster 2 investigation: cluster2_split_investigation.png\n")
cat("Cluster 2 temporal pattern: cluster2_temporal_pattern.png\n\n")
cat("=== KEY INSIGHTS ===\n")
cat("Check the diagnostic output above to understand:\n")
cat("- Which variables have the strongest unique effects\n")
cat("- Whether variables are independent or redundant\n")
cat("- How species (phytoplankton groups) align with environmental gradients\n")
cat("- Why Cluster 2 (diatoms) splits into different environmental conditions\n")