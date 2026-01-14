#==============================================================================
# CLEAN db-RDA WITH PROPER VARIABLE SELECTION
#==============================================================================

library(tidyverse)
library(here)
library(vegan)

#==============================================================================
# 1. LOAD DATA
#==============================================================================

chem <- read_csv(here("outputs", "chem_clusters_2025-12-03.csv"))
env <- read_csv(here("outputs", "enviro_2025-12-03.csv"))
# vert <- read_csv(here("outputs", "vertical_metrics.csv"))
vert <- read_csv(here("outputs", "vertical_structure_metrics.csv"))

#==============================================================================
# 2. DATA PREPARATION
#==============================================================================

chem_clean <- chem %>% select(-contains("..."))
env_clean <- env %>% select(-contains("..."))

vert <- vert %>% select(date, int_chla, vertical_concentration_index,
                        scm_intensity)



env_averaged <- env_clean %>%
  group_by(date) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  left_join(vert, by = "date")

rda_data <- chem_clean %>%
  left_join(env_averaged, by = "date") %>%
  drop_na()

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
    
    # Vertical gradients
    # delta_temp_5_30,
    delta_sigma_2_30,
    
    # Vertical structure
    int_chla,
    # centroid_depth,
    # profile_shape,
    # surface_depletion,
    # scm_depth,
    # vertical_gradient,
    vertical_concentration_index,
    scm_intensity,
    
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
# 6. PLOT
#==============================================================================
cluster_colors <- c("1" = "#DDAA33", "2" = "#228833", "3" = "#CC6677", "4" = "#4477AA")
cluster_labels <- c("1" = "1 (Mixed)", "2" = "2 (Diatoms)", 
                    "3" = "3 (Dinoflagellates)", "4" = "4 (Winter)")
plot_data <- as.data.frame(site_scores) %>%
  mutate(cluster = factor(rda_data$abs_clust))

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
    # Use absolute positioning only for scm_intensity
    label_x = if_else(variable == "scm_intensity", 0.2, CAP1 * 2.9),
    label_y = if_else(variable == "scm_intensity", -0.7, CAP2 * 2.9)
  )

var_exp <- summary(dbrda_final)$cont$importance[2, 1:2] * 100

p_dbrda <- ggplot() +
  geom_point(data = plot_data,
             aes(x = CAP1, y = CAP2, fill = cluster, shape = cluster),
             size = 7, alpha = 0.85, color = "black", stroke = 1) +
  
  geom_segment(data = env_arrows,
               aes(x = 0, y = 0, xend = CAP1 * 2.5, yend = CAP2 * 2.5),
               arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
               color = "black", linewidth = 1.8) +
  
  geom_text(data = env_arrows,
            aes(x = label_x, y = label_y, label = label),
            color = "black", size = 6, fontface = "bold") +
  
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
                     " | ", length(sig_vars), " significant variables")
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
ggsave(here("figures", "dbrda_final.png"), width = 14, height = 10, dpi = 400, bg = "white")

cat("\n✅ Analysis complete!\n")
cat("Final model has", length(sig_vars), "significant variables (p < 0.05)\n")
