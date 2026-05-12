#=============================================================================
# CORRELATION ANALYSIS: FIVE GENERA WITH WEIGHT
#=============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("WEIGHT CORRELATION ANALYSIS FOR FIVE GENERA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Define the five genera
five_genera <- c("Lactococcus", "Duncaniella", "Taurinivorans", "Lepagella", "CAG-873")

# Colors
diet_colors <- c("ND" = "#2E86AB", "HFD" = "#A23B72")
diet_fills <- c("ND" = "#2E86AB40", "HFD" = "#A23B7240")

# ============================================================================
# EXTRACT AND PREPARE DATA
# ============================================================================

# Build combined dataframe for five genera with weight
weight_cor_data <- data.frame()

for(g in five_genera) {
    abund <- get_genus_data(g)
    if(!is.null(abund)) {
        temp_df <- data.frame(
            Sample = names(abund),
            Genus = g,
            Abundance = as.numeric(abund),
            Day = metadata_clean$Day,
            Diet = metadata_clean$group1,
            Phase = metadata_clean$Phase,
            Cage = metadata_clean$cage,
            Weight = metadata_clean$weight,
            stringsAsFactors = FALSE
        )
        weight_cor_data <- rbind(weight_cor_data, temp_df)
    }
}

cat("Data dimensions:", nrow(weight_cor_data), "observations for", 
    length(unique(weight_cor_data$Genus)), "genera\n\n")

# ============================================================================
# CALCULATE CORRELATIONS: OVERALL, BY DIET, AND PARTIAL
# ============================================================================

# Initialize results dataframe
cor_results <- data.frame()

for(g in five_genera) {
    g_data <- weight_cor_data %>% filter(Genus == g)
    
    if(nrow(g_data) < 5) next
    
    cat("Processing:", g, "\n")
    cat("  Observations:", nrow(g_data), "(ND:", sum(g_data$Diet == "ND"), 
        ", HFD:", sum(g_data$Diet == "HFD"), ")\n")
    
    # ========================================================================
    # 1. OVERALL CORRELATION (all samples)
    # ========================================================================
    
    # Spearman correlation
    overall_spearman <- cor.test(g_data$Abundance, g_data$Weight, 
                                 method = "spearman", exact = FALSE)
    
    # Pearson correlation on log-transformed abundance
    overall_pearson <- cor.test(log10(g_data$Abundance + 1e-6), g_data$Weight, 
                                method = "pearson")
    
    cat("  Overall Spearman: ρ =", round(overall_spearman$estimate, 3), 
        ", p =", format(overall_spearman$p.value, scientific = TRUE, digits = 3), "\n")
    cat("  Overall Pearson: r =", round(overall_pearson$estimate, 3), 
        ", p =", format(overall_pearson$p.value, scientific = TRUE, digits = 3), "\n")
    
    # ========================================================================
    # 2. PER-DIET CORRELATION
    # ========================================================================
    
    for(diet in c("ND", "HFD")) {
        diet_data <- g_data %>% filter(Diet == diet)
        
        if(nrow(diet_data) >= 5) {
            # Spearman
            diet_spearman <- cor.test(diet_data$Abundance, diet_data$Weight, 
                                      method = "spearman", exact = FALSE)
            
            # Pearson on log-transformed
            diet_pearson <- cor.test(log10(diet_data$Abundance + 1e-6), diet_data$Weight, 
                                     method = "pearson")
            
            cor_results <- rbind(cor_results, data.frame(
                Genus = g,
                Diet = diet,
                N = nrow(diet_data),
                Spearman_rho = round(diet_spearman$estimate, 4),
                Spearman_p = diet_spearman$p.value,
                Pearson_r = round(diet_pearson$estimate, 4),
                Pearson_p = diet_pearson$p.value,
                Mean_Abundance = mean(diet_data$Abundance),
                Mean_Weight = mean(diet_data$Weight),
                stringsAsFactors = FALSE
            ))
            
            cat("  ", diet, "- Spearman: ρ =", round(diet_spearman$estimate, 3),
                ", p =", format(diet_spearman$p.value, scientific = TRUE, digits = 3), "\n")
        }
    }
    
    # ========================================================================
    # 3. PARTIAL CORRELATION (controlling for Day)
    # ========================================================================
    
    if(nrow(g_data) >= 10) {
        tryCatch({
            # Residuals from regressing Abundance on Day
            resid_abund <- residuals(lm(Abundance ~ Day, data = g_data))
            # Residuals from regressing Weight on Day
            resid_weight <- residuals(lm(Weight ~ Day, data = g_data))
            
            partial_cor <- cor.test(resid_abund, resid_weight, method = "spearman")
            
            cat("  Partial (Day controlled) - Spearman: ρ =", round(partial_cor$estimate, 3),
                ", p =", format(partial_cor$p.value, scientific = TRUE, digits = 3), "\n")
        }, error = function(e) {
            cat("  Partial correlation failed:", e$message, "\n")
        })
    }
    
    cat("\n")
}

# ============================================================================
# ADJUST P-VALUES
# ============================================================================

# Multiple testing correction (Benjamini-Hochberg)
cor_results <- cor_results %>%
    group_by(Diet) %>%
    mutate(
        Spearman_padj = p.adjust(Spearman_p, method = "BH"),
        Pearson_padj = p.adjust(Pearson_p, method = "BH")
    ) %>%
    ungroup() %>%
    mutate(
        # Add significance stars
        Significance = case_when(
            Spearman_padj < 0.001 ~ "***",
            Spearman_padj < 0.01  ~ "**",
            Spearman_padj < 0.05  ~ "*",
            TRUE ~ "ns"
        ),
        # Create label for plots
        Label = paste0(
            "ρ = ", round(Spearman_rho, 3),
            ifelse(Spearman_padj < 0.05, 
                   paste0("\np.adj = ", format(Spearman_padj, scientific = TRUE, digits = 2)),
                   "\nns")
        )
    )

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("FINAL CORRELATION RESULTS (WITH P-ADJUSTMENT)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

print(cor_results[, c("Genus", "Diet", "N", "Spearman_rho", "Spearman_padj", 
                       "Pearson_r", "Pearson_padj", "Significance")])

write.csv(cor_results, "five_genera_weight_correlations.csv", row.names = FALSE)

# ============================================================================
# OVERALL CORRELATION RESULTS (COMBINING BOTH DIETS)
# ============================================================================

overall_cor_results <- data.frame()

for(g in five_genera) {
    g_data <- weight_cor_data %>% filter(Genus == g)
    
    if(nrow(g_data) < 5) next
    
    overall_spearman <- cor.test(g_data$Abundance, g_data$Weight, 
                                 method = "spearman", exact = FALSE)
    overall_pearson <- cor.test(log10(g_data$Abundance + 1e-6), g_data$Weight, 
                                method = "pearson")
    
    overall_cor_results <- rbind(overall_cor_results, data.frame(
        Genus = g,
        N_total = nrow(g_data),
        N_ND = sum(g_data$Diet == "ND"),
        N_HFD = sum(g_data$Diet == "HFD"),
        Spearman_rho = round(overall_spearman$estimate, 4),
        Spearman_p = overall_spearman$p.value,
        Pearson_r = round(overall_pearson$estimate, 4),
        Pearson_p = overall_pearson$p.value,
        stringsAsFactors = FALSE
    ))
}

overall_cor_results <- overall_cor_results %>%
    mutate(
        Spearman_padj = p.adjust(Spearman_p, method = "BH"),
        Pearson_padj = p.adjust(Pearson_p, method = "BH"),
        Significance = case_when(
            Spearman_padj < 0.001 ~ "***",
            Spearman_padj < 0.01  ~ "**",
            Spearman_padj < 0.05  ~ "*",
            TRUE ~ "ns"
        )
    ) %>%
    arrange(Spearman_padj)

cat("\nOVERALL CORRELATIONS (Both diets combined):\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
print(overall_cor_results[, c("Genus", "N_total", "Spearman_rho", "Spearman_padj", 
                               "Pearson_r", "Pearson_padj", "Significance")])

write.csv(overall_cor_results, "five_genera_weight_correlations_overall.csv", row.names = FALSE)

# ============================================================================
# PLOT 1: SCATTER PLOTS WITH CORRELATION (FACETED BY GENUS)
# ============================================================================

p_scatter_facet <- ggplot(weight_cor_data, aes(x = Abundance, y = Weight)) +
    geom_point(aes(color = Diet), alpha = 0.6, size = 2) +
    geom_smooth(aes(color = Diet, fill = Diet), method = "lm", se = TRUE, alpha = 0.2) +
    facet_wrap(~ Genus, scales = "free_x", ncol = 3) +
    scale_color_manual(values = diet_colors, name = "Diet") +
    scale_fill_manual(values = diet_colors, name = "Diet") +
    scale_x_log10(labels = scales::scientific) +
    labs(
        title = "Weight Association of Differential Genera",
        subtitle = "Linear regression with 95% CI | X-axis: log10 scale",
        x = "Relative Abundance (log10 scale)",
        y = "Weight (g)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        strip.text = element_text(face = "italic", size = 12),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
    ) +
    # Add correlation labels
    geom_text(
        data = cor_results,
        aes(x = Inf, y = Inf, label = paste0(
            Diet, ": ", Label
        )),
        hjust = 1.05, vjust = 2,
        size = 3,
        inherit.aes = FALSE
    )

ggsave("weight_correlation_scatter_facet.pdf", p_scatter_facet, width = 14, height = 10)
print(p_scatter_facet)

# ============================================================================
# PLOT 2: SCATTER PLOTS BY DIET (COLORED BY GENUS)
# ============================================================================

genus_colors <- c("Lactococcus" = "#E41A1C", 
                  "Duncaniella" = "#377EB8", 
                  "Taurinivorans" = "#4DAF4A", 
                  "Lepagella" = "#984EA3", 
                  "CAG-873" = "#FF7F00")

p_scatter_diet <- ggplot(weight_cor_data, aes(x = Abundance, y = Weight, color = Genus)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(aes(fill = Genus), method = "lm", se = TRUE, alpha = 0.1) +
    facet_wrap(~ Diet, ncol = 2) +
    scale_color_manual(values = genus_colors, name = "Genus") +
    scale_fill_manual(values = genus_colors, guide = "none") +
    scale_x_log10(labels = scales::scientific) +
    labs(
        title = "Weight Association by Diet",
        subtitle = "Each genus fitted separately",
        x = "Relative Abundance (log10 scale)",
        y = "Weight (g)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(face = "italic"),
        plot.title = element_text(face = "bold")
    )

ggsave("weight_correlation_by_diet.pdf", p_scatter_diet, width = 12, height = 6)
print(p_scatter_diet)

# ============================================================================
# PLOT 3: CORRELATION COEFFICIENT COMPARISON
# ============================================================================

# Prepare data for coefficient plot
coef_plot_data <- cor_results %>%
    mutate(
        Genus_Diet = paste(Genus, Diet, sep = " - "),
        Direction = ifelse(Spearman_rho > 0, "Positive", "Negative")
    )

p_coef <- ggplot(coef_plot_data, aes(x = Spearman_rho, y = reorder(Genus_Diet, Spearman_rho))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_point(aes(color = Direction, size = -log10(Spearman_padj + 1e-300)), alpha = 0.8) +
    geom_errorbarh(aes(xmin = Spearman_rho - 0.05, xmax = Spearman_rho + 0.05), 
                   height = 0.2, linewidth = 0.8) +
    facet_wrap(~ Diet, scales = "free_y", ncol = 1) +
    scale_color_manual(
        values = c("Positive" = "#A23B72", "Negative" = "#2E86AB"),
        name = "Association"
    ) +
    scale_size_continuous(name = "-log10(p.adj)") +
    labs(
        title = "Spearman Correlation with Weight",
        subtitle = "Error bars: ±0.05 | Size: significance level",
        x = "Spearman ρ",
        y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.y = element_text(face = "italic", size = 10),
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
    ) +
    # Add correlation values
    geom_text(
        aes(label = paste0("ρ = ", round(Spearman_rho, 3), " ", Significance)),
        hjust = ifelse(coef_plot_data$Spearman_rho > 0, -0.3, 1.3),
        size = 3.5
    )

ggsave("weight_correlation_coefficients.pdf", p_coef, width = 10, height = 8)
print(p_coef)

# ============================================================================
# PLOT 4: OVERALL CORRELATION HEATMAP
# ============================================================================

# Create matrix for heatmap
cor_matrix_data <- overall_cor_results %>%
    select(Genus, Spearman_rho, Spearman_padj, Pearson_r, Pearson_padj)

# Prepare annotation
cor_annotation <- overall_cor_results %>%
    mutate(
        Direction = ifelse(Spearman_rho > 0, "Positive", "Negative"),
        Sig = case_when(
            Spearman_padj < 0.001 ~ "***",
            Spearman_padj < 0.01 ~ "**",
            Spearman_padj < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    )

p_heatmap_cor <- ggplot(cor_annotation, aes(x = "Weight", y = reorder(Genus, Spearman_rho))) +
    geom_tile(aes(fill = Spearman_rho), color = "white", size = 1) +
    geom_text(aes(label = paste0(round(Spearman_rho, 3), "\n", Sig)),
              size = 4, fontface = "bold") +
    scale_fill_gradient2(
        low = "#2E86AB", mid = "white", high = "#A23B72",
        midpoint = 0, limits = c(-1, 1),
        name = "Spearman ρ"
    ) +
    labs(
        title = "Weight Correlation Summary",
        subtitle = "Overall (both diets combined)",
        x = "",
        y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.y = element_text(face = "italic", size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold"),
        panel.grid = element_blank()
    )

ggsave("weight_correlation_heatmap.pdf", p_heatmap_cor, width = 6, height = 5)
print(p_heatmap_cor)

# ============================================================================
# PLOT 5: COMBINED SCATTER WITH HISTOGRAMS
# ============================================================================

plot_list <- list()

for(g in five_genera) {
    g_data <- weight_cor_data %>% filter(Genus == g)
    if(nrow(g_data) < 5) next
    
    cor_info <- overall_cor_results %>% filter(Genus == g)
    
    # Main scatter plot
    p_main <- ggplot(g_data, aes(x = Abundance, y = Weight, color = Diet)) +
        geom_point(alpha = 0.7, size = 2.5) +
        geom_smooth(aes(fill = Diet), method = "lm", se = TRUE, alpha = 0.2) +
        scale_color_manual(values = diet_colors) +
        scale_fill_manual(values = diet_colors) +
        scale_x_log10(labels = scales::scientific) +
        labs(
            x = "Relative Abundance",
            y = "Weight (g)",
            title = bquote(italic(.(g)))
        ) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none") +
        annotate("text", x = Inf, y = Inf,
                 label = paste0("ρ = ", round(cor_info$Spearman_rho, 3),
                                " ", cor_info$Significance,
                                "\np.adj = ", format(cor_info$Spearman_padj, 
                                                     scientific = TRUE, digits = 2)),
                 hjust = 1.1, vjust = 1.5, size = 3.5)
    
    # Marginal density plot for abundance
    p_dens_x <- ggplot(g_data, aes(x = Abundance, fill = Diet)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = diet_colors) +
        scale_x_log10(labels = scales::scientific) +
        theme_void() +
        theme(legend.position = "none")
    
    # Marginal density plot for weight
    p_dens_y <- ggplot(g_data, aes(x = Weight, fill = Diet)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = diet_colors) +
        theme_void() +
        theme(legend.position = "none") +
        coord_flip()
    
    # Combine using patchwork
    combined <- p_dens_x + plot_spacer() + p_main + p_dens_y +
        plot_layout(
            widths = c(1, 0.1),
            heights = c(0.2, 1)
        )
    
    plot_list[[g]] <- combined
}

if(length(plot_list) > 0) {
    all_scatter_plots <- wrap_plots(plot_list, ncol = 3)
    ggsave("weight_correlation_detailed_scatter.pdf", all_scatter_plots, 
           width = 16, height = 8)
}

# ============================================================================
# SUMMARY TABLE AND REPORT
# ============================================================================

sink("weight_correlation_report.txt")

cat("========================================\n")
cat("WEIGHT CORRELATION ANALYSIS REPORT\n")
cat("Five High-Confidence Differential Genera\n")
cat("========================================\n\n")
cat("Date:", date(), "\n\n")

cat("OVERALL CORRELATIONS (ND + HFD combined):\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
for(i in 1:nrow(overall_cor_results)) {
    cat(sprintf("%-15s: ρ = %+7.4f | p.adj = %8s | %s\n",
                overall_cor_results$Genus[i],
                overall_cor_results$Spearman_rho[i],
                format(overall_cor_results$Spearman_padj[i], scientific = TRUE, digits = 3),
                overall_cor_results$Significance[i]))
}

cat("\nPER-DIET CORRELATIONS:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
for(i in 1:nrow(cor_results)) {
    cat(sprintf("%-15s [%3s]: ρ = %+7.4f | p.adj = %8s | %s | N = %d\n",
                cor_results$Genus[i],
                cor_results$Diet[i],
                cor_results$Spearman_rho[i],
                format(cor_results$Spearman_padj[i], scientific = TRUE, digits = 3),
                cor_results$Significance[i],
                cor_results$N[i]))
}

cat("\nKEY FINDINGS:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

sig_overall <- overall_cor_results %>% filter(Spearman_padj < 0.05)
if(nrow(sig_overall) > 0) {
    cat("Significantly correlated with weight (overall):\n")
    for(i in 1:nrow(sig_overall)) {
        direction <- ifelse(sig_overall$Spearman_rho[i] > 0, "positively", "negatively")
        cat("  -", sig_overall$Genus[i], ":", direction, "correlated\n")
    }
} else {
    cat("No genera significantly correlated with weight after p-value adjustment.\n")
}

sig_by_diet <- cor_results %>% filter(Spearman_padj < 0.05)
if(nrow(sig_by_diet) > 0) {
    cat("\nSignificant diet-specific correlations:\n")
    for(i in 1:nrow(sig_by_diet)) {
        direction <- ifelse(sig_by_diet$Spearman_rho[i] > 0, "positively", "negatively")
        cat("  -", sig_by_diet$Genus[i], "(", sig_by_diet$Diet[i], "):", direction, "correlated\n")
    }
}

sink()

# ============================================================================
# SAVE RESULTS
# ============================================================================

save(weight_cor_data, cor_results, overall_cor_results,
     file = "weight_correlation_analysis.RData")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("WEIGHT CORRELATION ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
cat("Files saved:\n")
cat("  - five_genera_weight_correlations.csv (per-diet results)\n")
cat("  - five_genera_weight_correlations_overall.csv (overall results)\n")
cat("  - weight_correlation_scatter_facet.pdf\n")
cat("  - weight_correlation_by_diet.pdf\n")
cat("  - weight_correlation_coefficients.pdf\n")
cat("  - weight_correlation_heatmap.pdf\n")
cat("  - weight_correlation_detailed_scatter.pdf\n")
cat("  - weight_correlation_report.txt\n")
cat("  - weight_correlation_analysis.RData\n")

# Print quick summary
cat("\nQuick Summary:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for(i in 1:nrow(overall_cor_results)) {
    cat(sprintf("  %-15s: ρ = %+7.4f [%s]\n",
                overall_cor_results$Genus[i],
                overall_cor_results$Spearman_rho[i],
                overall_cor_results$Significance[i]))
}