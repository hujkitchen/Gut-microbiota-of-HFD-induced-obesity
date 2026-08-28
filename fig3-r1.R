#=============================================================================
# FIGURE 3 (REVISED): PER-MOUSE AGGREGATED CORRELATIONS
# Correlation between genus abundance and body weight using per-mouse means
#=============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggrepel)
library(patchwork)

cat("Figure 3: Per-Mouse Aggregated Correlations\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

#=============================================================================
# 1. PREPARE PER-MOUSE DATA
#=============================================================================

# Define the five genera
five_genera <- c("Lactococcus", "Duncaniella", "Taurinivorans", "Lepagella", "CAG-873")


if("mice" %in% colnames(metadata_clean)) {
    mouse_id <- metadata_clean$mice
} else {
    mouse_id <- gsub("_DAY.*", "", metadata_clean$SampleID)
}

five_genera_clean <- gsub("^g__", "", five_genera)

per_mouse_data <- data.frame()

for(g in five_genera_clean) {
    if(g %in% rownames(genus_rel)) {
        abund <- as.numeric(genus_rel[g, ])
    } else if(paste0("g__", g) %in% rownames(genus_rel)) {
        abund <- as.numeric(genus_rel[paste0("g__", g), ])
    } else {
        next
    }
    
    mouse_df <- data.frame(
        Mouse = factor(mouse_id),
        Diet = metadata_clean$group1,
        Cage = metadata_clean$cage,
        Abundance = abund,
        Weight = as.numeric(metadata_clean$weight)
    )
    
    mouse_agg <- mouse_df %>%
        group_by(Mouse, Diet, Cage) %>%
        summarise(
            MeanAbundance = mean(Abundance, na.rm = TRUE),
            MeanWeight = mean(Weight, na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        mutate(Genus = g)
    
    per_mouse_data <- rbind(per_mouse_data, mouse_agg)
}

cat("Per-mouse data:", nrow(per_mouse_data), "rows\n")

#=============================================================================
# 2. CALCULATE CORRELATIONS
#=============================================================================

cor_results <- data.frame()

for(g in five_genera_clean) {
    g_data <- per_mouse_data %>% filter(Genus == g)
    
    cor_overall <- cor.test(g_data$MeanAbundance, g_data$MeanWeight,
                            method = "spearman", exact = FALSE)
    
    cor_results <- rbind(cor_results, data.frame(
        Genus = g,
        Overall_rho = round(cor_overall$estimate, 3),
        Overall_p = cor_overall$p.value,
        stringsAsFactors = FALSE
    ))
}

# Adjust p-values
cor_results$Overall_padj <- p.adjust(cor_results$Overall_p, method = "BH")

# Add significance and direction
cor_results <- cor_results %>%
    mutate(
        Sig = case_when(
            Overall_padj < 0.001 ~ "***",
            Overall_padj < 0.01 ~ "**",
            Overall_padj < 0.05 ~ "*",
            TRUE ~ "ns"
        ),
        Direction = ifelse(Overall_rho > 0, "Positive", "Negative")
    )

cat("\nCorrelation results:\n")
print(cor_results)

#=============================================================================
# 3. SCATTER PLOTS (Panel A)
#=============================================================================

diet_colors <- c("ND" = "#2E86AB", "HFD" = "#A23B72")

scatter_plots <- list()

for(g in five_genera_clean) {
    g_data <- per_mouse_data %>% filter(Genus == g)
    stats <- cor_results %>% filter(Genus == g)
    
    anno_text <- paste0("ρ = ", stats$Overall_rho, " ", stats$Sig)
    
    scatter_plots[[g]] <- ggplot(g_data, aes(x = MeanAbundance, y = MeanWeight)) +
        geom_point(aes(color = Diet), size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "grey30", 
                    linewidth = 0.8, alpha = 0.15) +
        geom_smooth(aes(color = Diet), method = "lm", se = TRUE,
                    linewidth = 0.6, alpha = 0.15) +
        scale_color_manual(values = diet_colors) +
        scale_x_log10(labels = scales::scientific) +
        annotate("text", x = Inf, y = Inf, label = anno_text,
                 hjust = 1.1, vjust = 1.1, size = 4, fontface = "bold") +
        labs(title = bquote(italic(.(g))),
             x = "Mean Relative Abundance",
             y = "Mean Body Weight (g)") +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(face = "bold.italic", size = 13),
              legend.position = "none")
}

p_scatter <- wrap_plots(scatter_plots, ncol = 3) +
    plot_annotation(
        title = "A. Per-Mouse Scatter Plots",
        subtitle = paste0("n = ", length(unique(mouse_id)), " mice | Each point = one mouse"),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey40")
        )
    )

#=============================================================================
# 4. FOREST PLOT (Panel B) - FIXED
#=============================================================================

# Ensure Direction exists and is character
forest_data <- cor_results %>%
    mutate(
        Genus = factor(Genus, levels = Genus[order(Overall_rho)]),
        Direction = as.character(Direction),
        Label = paste0("ρ = ", Overall_rho, " ", Sig)
    )

cat("\nForest data columns:", paste(colnames(forest_data), collapse = ", "), "\n")
cat("Direction values:", paste(forest_data$Direction, collapse = ", "), "\n")

p_forest <- ggplot(forest_data, aes(x = Overall_rho, y = Genus)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = c(-1, -0.5, 0.5, 1), linetype = "dotted", 
               color = "grey85") +
    # Points with explicit color mapping
    geom_point(aes(color = Direction), size = 6, alpha = 0.9) +
    # Labels
    geom_text(aes(label = Label), hjust = -0.3, size = 4, fontface = "bold") +
    # Color scale
    scale_color_manual(
        values = c("Positive" = "#A23B72", "Negative" = "#2E86AB"),
        guide = "none"
    ) +
    labs(
        title = "B. Correlation Coefficients",
        subtitle = "Spearman ρ | BH-corrected p-values",
        x = "Spearman ρ",
        y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.y = element_text(face = "italic", size = 13),
        plot.title = element_text(face = "bold", size = 13)
    ) +
    xlim(-1.3, 1.3)

print(p_forest)

#=============================================================================
# 5. COMBINED FIGURE
#=============================================================================

p_combined <- p_scatter / p_forest +
    plot_annotation(
        title = "Figure 3. Weight Associations of Five Focal Genera (Per-Mouse Analysis)",
        tag_levels = "A",
        theme = theme(plot.title = element_text(face = "bold", size = 15))
    ) +
    plot_layout(heights = c(1.5, 1))

ggsave("Figure3_per_mouse_combined.pdf", p_combined,
       width = 14, height = 12, dpi = 300)
print(p_combined)

#=============================================================================
# 6. SAVE
#=============================================================================

write.csv(cor_results, "Table_per_mouse_correlations.csv", row.names = FALSE)
save(per_mouse_data, cor_results, p_scatter, p_forest, p_combined,
     file = "figure3_per_mouse_results.RData")

cat("\nDone!\n")