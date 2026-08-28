#=============================================================================
# R2-4: REPLACE VENN DIAGRAM WITH EFFECT SIZE SCATTER
#=============================================================================


library(ggplot2)
library(dplyr)
library(ggrepel)

cat("Creating Figure 2A: Method Concordance Plot (Fixed)\n\n")

#=============================================================================
# 1. PREPARE DATA - Handle Duplicates
#=============================================================================

# Clean five_genera
five_genera_clean <- gsub("^g__", "", five_genera)

# Inspect MaAsLin3 data structure to understand duplicates
cat("MaAsLin3 columns:", paste(colnames(maaslin_df), collapse = ", "), "\n")
cat("MaAsLin3 unique names:", length(unique(maaslin_df$name)), "\n")
cat("MaAsLin3 metadata values:", paste(unique(maaslin_df$metadata), collapse = ", "), "\n")

# Check if there are multiple model types (abundance vs prevalence)
if("model" %in% colnames(maaslin_df)) {
    cat("\nModel types in MaAsLin3 results:\n")
    print(table(maaslin_df$model))
} else if("name" %in% colnames(maaslin_df)) {
    cat("\nName values in MaAsLin3 results:\n")
    print(table(maaslin_df$name))
}

#=============================================================================
# 2. SELECT BEST RESULT PER GENUS
#=============================================================================

# For MaAsLin3: keep the most significant result per genus (lowest q-value)
# or filter for a specific model type if available

if("model" %in% colnames(maaslin_df)) {
    cat("\nFiltering MaAsLin3 results by model type...\n")
    
    # Prefer abundance model over prevalence
    if("abundance" %in% maaslin_df$model) {
        maaslin_filtered <- maaslin_df %>%
            filter(model == "abundance")
        cat("Using abundance model results\n")
    } else {
        maaslin_filtered <- maaslin_df
        cat("Using all models (no abundance model found)\n")
    }
} else {
    maaslin_filtered <- maaslin_df
}

# For each genus, keep the row with the lowest q-value
maaslin_unique <- maaslin_filtered %>%
    group_by(Genus_clean) %>%
    slice_min(order_by = MaAsLin_q, n = 1, with_ties = FALSE) %>%
    ungroup()

cat("\nMaAsLin3 unique genera:", nrow(maaslin_unique), "\n")

# ANCOM-BC2: keep unique per genus
ancom_unique <- ancom_df %>%
    group_by(Genus_clean) %>%
    slice_min(order_by = ANCOM_q, n = 1, with_ties = FALSE) %>%
    ungroup()

cat("ANCOM-BC2 unique genera:", nrow(ancom_unique), "\n")

#=============================================================================
# 3. MERGE DATA PROPERLY
#=============================================================================

comparison_df <- ancom_unique %>%
    full_join(maaslin_unique, by = "Genus_clean") %>%
    mutate(
        Detected_By = case_when(
            !is.na(ANCOM_LFC) & !is.na(MaAsLin_coef) ~ "Both methods",
            !is.na(ANCOM_LFC) & is.na(MaAsLin_coef) ~ "ANCOM-BC2 only",
            is.na(ANCOM_LFC) & !is.na(MaAsLin_coef) ~ "MaAsLin3 only",
            TRUE ~ "Neither"
        ),
        is_focal = Genus_clean %in% five_genera_clean,
        Label = ifelse(is_focal, Genus_clean, "")
    )

cat("\nDetection summary (unique genera):\n")
print(table(comparison_df$Detected_By))

cat("\nFocal genera (unique):\n")
print(comparison_df %>% filter(is_focal))

#=============================================================================
# 4. CORRECT ANCOM-BC2 DIRECTION
#=============================================================================

# ANCOM-BC2 column was lfc_group1ND = log(ND/HFD)
# Need to negate for HFD vs ND
comparison_df <- comparison_df %>%
    mutate(
        ANCOM_LFC_corrected = -ANCOM_LFC  # Negate for HFD vs ND
    )

cat("\nANCOM-BC2 LFC corrected (negated for HFD vs ND):\n")
print(comparison_df %>% filter(is_focal) %>% 
        select(Genus_clean, ANCOM_LFC, ANCOM_LFC_corrected, MaAsLin_coef))

#=============================================================================
# 5. CREATE MAIN SCATTER PLOT
#=============================================================================

detection_colors <- c(
    "Both methods" = "#E63946",
    "ANCOM-BC2 only" = "#2E86AB",
    "MaAsLin3 only" = "#4CAF50",
    "Neither" = "grey80"
)

# Calculate concordance correlation
concordant_data <- comparison_df %>% filter(Detected_By == "Both methods")
if(nrow(concordant_data) >= 3) {
    cor_test <- cor.test(concordant_data$ANCOM_LFC_corrected, 
                         concordant_data$MaAsLin_coef, 
                         method = "spearman")
    cat("\nSpearman concordance (corrected): rho =", round(cor_test$estimate, 3),
        ", p =", format(cor_test$p.value, scientific = TRUE, digits = 2), "\n")
}

p_concordance <- ggplot(comparison_df, aes(x = ANCOM_LFC_corrected, y = MaAsLin_coef)) +
    # Reference lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                color = "darkgreen", linewidth = 0.7, alpha = 0.5) +
    
    # Points
    geom_point(
        aes(color = Detected_By, 
            size = Detected_By == "Both methods",
            alpha = Detected_By == "Both methods"),
        stroke = 0.3
    ) +
    
    # Labels for focal genera
    geom_text_repel(
        data = comparison_df %>% filter(is_focal),
        aes(label = Genus_clean),
        fontface = "bold.italic",
        size = 5,
        color = "#8B0000",
        box.padding = 0.8,
        point.padding = 0.5,
        max.overlaps = 10,
        segment.color = "#8B0000",
        segment.size = 0.5
    ) +
    
    # Scales
    scale_color_manual(values = detection_colors, name = "Detected by", drop = FALSE) +
    scale_size_manual(values = c("TRUE" = 4.5, "FALSE" = 2), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.45), guide = "none") +
    
    # Labels
    labs(
        title = "Method Concordance: ANCOM-BC2 vs MaAsLin3",
        subtitle = paste0(
            "ANCOM-BC2 LFC corrected (negated for HFD vs ND)\n",
            "Red = detected by both methods | Dashed line = perfect agreement"
        ),
        x = "ANCOM-BC2 Log2 Fold Change (HFD vs ND)",
        y = "MaAsLin3 Coefficient (HFD vs ND)",
        caption = paste0(
            "Both: ", sum(comparison_df$Detected_By == "Both methods"), " genera | ",
            "ANCOM-BC2 only: ", sum(comparison_df$Detected_By == "ANCOM-BC2 only"), " | ",
            "MaAsLin3 only: ", sum(comparison_df$Detected_By == "MaAsLin3 only")
        )
    ) +
    
    theme_bw(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10, color = "grey40", lineheight = 1.3),
        plot.caption = element_text(size = 9, color = "grey50"),
        legend.position = "bottom"
    )

ggsave("Revised_Figure2A_concordance.pdf", p_concordance, width = 10, height = 8, dpi = 300)
print(p_concordance)

#=============================================================================
# 6. UPSET-STYLE BAR PLOT
#=============================================================================

counts_df <- comparison_df %>%
    count(Detected_By) %>%
    filter(Detected_By != "Neither")

p_upset <- ggplot(counts_df, aes(x = Detected_By, y = n, fill = Detected_By)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.5) +
    geom_text(aes(label = n), vjust = -0.5, size = 6, fontface = "bold") +
    scale_fill_manual(values = detection_colors, guide = "none") +
    labs(
        title = "Differentially Abundant Genera by Detection Method",
        x = "",
        y = "Number of Genera"
    ) +
    theme_bw(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11)
    ) +
    ylim(0, max(counts_df$n) * 1.2)

ggsave("Revised_Figure2A_upset.pdf", p_upset, width = 6, height = 5, dpi = 300)
print(p_upset)

#=============================================================================
# 7. COMBINED FIGURE
#=============================================================================

p_combined <- p_concordance | p_upset +
    plot_annotation(
        title = "Differential Abundance: Method Concordance",
        tag_levels = "A",
        theme = theme(plot.title = element_text(face = "bold", size = 16))
    )

ggsave("Figure2A_combined.pdf", p_combined, width = 16, height = 8, dpi = 300)
print(p_combined)

#=============================================================================
# 8. SAVE
#=============================================================================

write.csv(comparison_df, "Table_method_concordance.csv", row.names = FALSE)
save(comparison_df, p_concordance, p_upset, p_combined,
     file = "figure2A_results.RData")

cat("\nDone!\n")