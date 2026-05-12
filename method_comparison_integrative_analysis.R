#=============================================================================
# COMPREHENSIVE METHOD COMPARISON AND FINAL INTEGRATIVE ANALYSIS
# ANCOM-BC2 vs MaAsLin3
#=============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggvenn)
library(corrplot)
library(patchwork)

# ============================================================================
# LOAD AND HARMONIZE RESULTS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("COMPREHENSIVE METHOD COMPARISON: ANCOM-BC2 vs MaAsLin3\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Load ANCOM-BC2 results
if(file.exists("ancombc2_final_results.csv")) {
    ancom_results <- read.csv("ancombc2_final_results.csv")
    cat("ANCOM-BC2 results loaded:", nrow(ancom_results), "genera\n")
} else if(exists("sig_ancom")) {
    ancom_results <- sig_ancom
    cat("Using ANCOM-BC2 results from workspace:", nrow(ancom_results), "genera\n")
}

# Load MaAsLin3 results
maaslin_file <- "maaslin3_diet_effect/significant_results.tsv"
if(file.exists(maaslin_file)) {
    maaslin_all <- read.delim(maaslin_file)
    maaslin_diet <- maaslin_all %>%
        filter(metadata == "group1") %>%
        mutate(
            taxon = gsub("^g__", "", feature),           # Remove g__ prefix
            taxon_clean = gsub("_", " ", taxon)           # Clean up names
        )
    cat("MaAsLin3 diet results loaded:", nrow(maaslin_diet), "genera\n")
}

# Clean ANCOM-BC2 names
ancom_results <- ancom_results %>%
    mutate(
        taxon_clean = gsub("^g__", "", taxon),
        taxon_clean = gsub("_", " ", taxon_clean)
    )

cat("\nANCOM-BC2 genera:", paste(ancom_results$taxon_clean, collapse = ", "), "\n")
cat("MaAsLin3 genera:", paste(maaslin_diet$taxon_clean, collapse = ", "), "\n")

# ============================================================================
# FIND OVERLAP WITH NAME HARMONIZATION
# ============================================================================

# Create sets with cleaned names
ancom_set <- unique(ancom_results$taxon_clean)
maaslin_set <- unique(maaslin_diet$taxon_clean)

# Find overlaps
overlap_genera <- intersect(ancom_set, maaslin_set)
ancom_only <- setdiff(ancom_set, maaslin_set)
maaslin_only <- setdiff(maaslin_set, ancom_set)

cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("METHOD COMPARISON RESULTS\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")
cat("Total ANCOM-BC2 significant genera:", length(ancom_set), "\n")
cat("Total MaAsLin3 significant genera:", length(maaslin_set), "\n")
cat("\nOverlap (detected by both):", length(overlap_genera), "\n")
cat("ANCOM-BC2 only:", length(ancom_only), "\n")
cat("MaAsLin3 only:", length(maaslin_only), "\n")

if(length(overlap_genera) > 0) {
    cat("\nGenera detected by BOTH methods:\n")
    cat(paste("  -", overlap_genera, "\n"))
}

if(length(ancom_only) > 0) {
    cat("\nGenera detected ONLY by ANCOM-BC2:\n")
    cat(paste("  -", ancom_only, "\n"))
}

if(length(maaslin_only) > 0) {
    cat("\nGenera detected ONLY by MaAsLin3:\n")
    cat(paste("  -", maaslin_only, "\n"))
}

# ============================================================================
# VENN DIAGRAM
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("VENN DIAGRAM\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Create Venn diagram data
venn_data <- list(
    `ANCOM-BC2` = ancom_set,
    `MaAsLin3` = maaslin_set
)

p_venn <- ggvenn(
    venn_data,
    fill_color = c("#2E86AB", "#A23B72"),
    stroke_size = 0.5,
    set_name_size = 5,
    text_size = 4
) +
    ggtitle("Overlap of Differential Genera: ANCOM-BC2 vs MaAsLin3")

ggsave("method_comparison_venn.pdf", p_venn, width = 8, height = 8)
print(p_venn)

# ============================================================================
# COMPARISON OF EFFECT SIZES
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("EFFECT SIZE COMPARISON\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Create comparison dataframe for overlapping genera
if(length(overlap_genera) > 0) {
    comparison_df <- data.frame(
        genus = overlap_genera,
        ancom_lfc = sapply(overlap_genera, function(g) {
            idx <- which(ancom_results$taxon_clean == g)
            if(length(idx) > 0) return(ancom_results$lfc[idx[1]])
            return(NA)
        }),
        maaslin_coef = sapply(overlap_genera, function(g) {
            idx <- which(maaslin_diet$taxon_clean == g)
            if(length(idx) > 0) return(maaslin_diet$coef[idx[1]])
            return(NA)
        }),
        ancom_qval = sapply(overlap_genera, function(g) {
            idx <- which(ancom_results$taxon_clean == g)
            if(length(idx) > 0) return(ancom_results$qval[idx[1]])
            return(NA)
        }),
        maaslin_qval = sapply(overlap_genera, function(g) {
            idx <- which(maaslin_diet$taxon_clean == g)
            if(length(idx) > 0) return(maaslin_diet$qval_individual[idx[1]])
            return(NA)
        })
    ) %>%
    mutate(
        concordant = sign(ancom_lfc) == sign(maaslin_coef),
        mean_effect = (ancom_lfc + maaslin_coef) / 2,
        mean_qval = -log10((ancom_qval + maaslin_qval) / 2)
    )
    
    cat("Effect size comparison for overlapping genera:\n")
    print(comparison_df[, c("genus", "ancom_lfc", "maaslin_coef", "concordant")])
    
    # Save comparison
    write.csv(comparison_df, "method_effect_comparison.csv", row.names = FALSE)
    
    # Concordance correlation
    if(nrow(comparison_df) >= 3) {
        concord_cor <- cor.test(comparison_df$ancom_lfc, comparison_df$maaslin_coef, 
                                method = "spearman")
        cat("\nSpearman correlation between effect sizes:", 
            round(concord_cor$estimate, 3), 
            "(p =", format(concord_cor$p.value, scientific = TRUE, digits = 2), ")\n")
        
        # Scatter plot comparing effect sizes
        p_concord <- ggplot(comparison_df, aes(x = ancom_lfc, y = maaslin_coef)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
            geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
            geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
            geom_point(aes(color = concordant, size = mean_qval), alpha = 0.8) +
            geom_smooth(method = "lm", se = TRUE, color = "darkred", alpha = 0.2) +
            geom_text_repel(
                data = comparison_df,
                aes(label = genus),
                size = 3.5,
                fontface = "italic",
                max.overlaps = 15,
                box.padding = 0.5
            ) +
            labs(
                x = "ANCOM-BC2 Log Fold Change",
                y = "MaAsLin3 Coefficient",
                title = "Effect Size Concordance Between Methods",
                subtitle = paste0("Spearman ρ = ", round(concord_cor$estimate, 3),
                                  ", p = ", format(concord_cor$p.value, scientific = TRUE, digits = 2))
            ) +
            scale_color_manual(
                values = c("TRUE" = "#4CAF50", "FALSE" = "#FF5722"),
                labels = c("Discordant", "Concordant"),
                name = "Direction"
            ) +
            scale_size_continuous(name = "Mean -log10(Q)") +
            theme_bw(base_size = 12) +
            theme(legend.position = "bottom")
        
        ggsave("method_effect_concordance.pdf", p_concord, width = 10, height = 8)
        print(p_concord)
    }
}

# ============================================================================
# CREATE HIGH-CONFIDENCE GENE SET
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("HIGH-CONFIDENCE GENERA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Genera detected by both methods (highest confidence)
high_confidence <- data.frame(
    genus = character(),
    evidence = character(),
    stringsAsFactors = FALSE
)

# Add concordant genera
if(length(overlap_genera) > 0) {
    high_confidence <- rbind(high_confidence, 
                             data.frame(genus = overlap_genera, 
                                        evidence = "Both methods",
                                        stringsAsFactors = FALSE))
}

# Add genera from single methods with trend information
if(length(ancom_only) > 0) {
    for(g in ancom_only) {
        idx <- which(ancom_results$taxon_clean == g)
        direction <- ifelse(ancom_results$lfc[idx[1]] > 0, 
                            "↑ HFD (ANCOM-BC2)", "↓ HFD (ANCOM-BC2)")
        high_confidence <- rbind(high_confidence,
                                 data.frame(genus = g,
                                            evidence = direction,
                                            stringsAsFactors = FALSE))
    }
}

if(length(maaslin_only) > 0) {
    for(g in maaslin_only) {
        idx <- which(maaslin_diet$taxon_clean == g)
        direction <- ifelse(maaslin_diet$coef[idx[1]] > 0,
                            "↑ HFD (MaAsLin3)", "↓ HFD (MaAsLin3)")
        high_confidence <- rbind(high_confidence,
                                 data.frame(genus = g,
                                            evidence = direction,
                                            stringsAsFactors = FALSE))
    }
}

# Add effect sizes
high_confidence$ancom_lfc <- sapply(high_confidence$genus, function(g) {
    idx <- which(ancom_results$taxon_clean == g)
    if(length(idx) > 0) return(ancom_results$lfc[idx[1]])
    return(NA)
})

high_confidence$maaslin_coef <- sapply(high_confidence$genus, function(g) {
    idx <- which(maaslin_diet$taxon_clean == g)
    if(length(idx) > 0) return(maaslin_diet$coef[idx[1]])
    return(NA)
})

cat("High-confidence differential genera:\n")
print(high_confidence)
write.csv(high_confidence, "high_confidence_genera.csv", row.names = FALSE)

# ============================================================================
# COMBINED HEATMAP OF ALL SIGNIFICANT GENERA
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("COMBINED HEATMAP\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Get all unique significant genera
all_sig_genera <- unique(c(ancom_set, maaslin_set))
cat("Total unique significant genera across methods:", length(all_sig_genera), "\n")

# Match with genus_rel
all_sig_in_data <- all_sig_genera[all_sig_genera %in% rownames(genus_rel)]

if(length(all_sig_in_data) > 0) {
    # Add original naming for lookup
    orig_names <- c(
        setNames(ancom_results$taxon, ancom_results$taxon_clean),
        setNames(maaslin_diet$feature, maaslin_diet$taxon_clean)
    )
    orig_names <- orig_names[!duplicated(names(orig_names))]
    
    # Get relative abundance data
    heatmap_data <- genus_rel[all_sig_in_data, , drop = FALSE]
    rownames(heatmap_data) <- all_sig_in_data  # Use clean names
    
    # Log transform
    heatmap_data_log <- log10(heatmap_data + 1e-5)
    
    # Create annotation
    anno_col <- data.frame(
        Diet = metadata_clean$group1,
        Phase = metadata_clean$Phase,
        row.names = colnames(heatmap_data_log)
    )
    
    anno_colors <- list(
        Diet = c("ND" = "#2E86AB", "HFD" = "#A23B72"),
        Phase = c("Acute" = "#440154", "Early" = "#3B528B", 
                  "Mid" = "#21908C", "Late" = "#5DC863")
    )
    
    # Create row annotation for method detection
    method_annotation <- data.frame(
        Detected_by = ifelse(all_sig_in_data %in% overlap_genera, "Both",
                             ifelse(all_sig_in_data %in% ancom_only, "ANCOM-BC2", "MaAsLin3")),
        row.names = all_sig_in_data
    )
    
    anno_colors$Detected_by <- c(
        "Both" = "#4CAF50",
        "ANCOM-BC2" = "#2E86AB",
        "MaAsLin3" = "#A23B72"
    )
    
    pdf("combined_significant_genera_heatmap.pdf", 
        width = 14, height = max(8, length(all_sig_in_data) * 0.4))
    
    if(length(all_sig_in_data) > 1) {
        pheatmap(heatmap_data_log,
                 scale = "row",
                 clustering_distance_rows = "correlation",
                 clustering_distance_cols = "euclidean",
                 annotation_col = anno_col,
                 annotation_row = method_annotation,
                 annotation_colors = anno_colors,
                 show_rownames = TRUE,
                 fontsize_row = 9,
                 main = paste("All Differential Genera (n =", length(all_sig_in_data), ")"),
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
    }
    dev.off()
    cat("Combined heatmap saved to combined_significant_genera_heatmap.pdf\n")
}

# ============================================================================
# TRAJECTORIES OF HIGH-CONFIDENCE GENERA
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("TRAJECTORIES OF HIGH-CONFIDENCE GENERA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Focus on genera detected by both methods or with strong effects
high_conf_genera <- high_confidence %>%
    filter(evidence == "Both methods" | grepl("↑|↓", evidence)) %>%
    pull(genus)

if(length(high_conf_genera) > 0) {
    traj_plots <- list()
    
    for(g in high_conf_genera) {
        if(g %in% rownames(genus_rel)) {
            traj_df <- data.frame(
                Abundance = genus_rel[g, ],
                Day = metadata_clean$Day,
                Diet = metadata_clean$group1,
                Phase = metadata_clean$Phase
            )
            
            traj_summary <- traj_df %>%
                group_by(Diet, Day) %>%
                summarise(
                    Mean = mean(Abundance),
                    SE = sd(Abundance) / sqrt(n()),
                    .groups = 'drop'
                )
            
            # Get statistics
            ancom_idx <- which(ancom_results$taxon_clean == g)
            maaslin_idx <- which(maaslin_diet$taxon_clean == g)
            
            stats_text <- paste0(
                ifelse(length(ancom_idx) > 0,
                       paste0("ANCOM LFC=", round(ancom_results$lfc[ancom_idx[1]], 2)), ""),
                ifelse(length(maaslin_idx) > 0,
                       paste0("\nMaAsLin coef=", round(maaslin_diet$coef[maaslin_idx[1]], 2)), "")
            )
            
            traj_plots[[g]] <- ggplot(traj_summary, aes(x = Day, y = Mean, color = Diet)) +
                geom_line(size = 1.2) +
                geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE, fill = Diet), 
                            alpha = 0.2, color = NA) +
                labs(
                    x = "Day", 
                    y = "Relative Abundance",
                    title = bquote(italic(.(g)))
                ) +
                annotate("text", x = max(traj_summary$Day) * 0.7, 
                         y = max(traj_summary$Mean + traj_summary$SE) * 0.9,
                         label = stats_text, size = 2.5, hjust = 0) +
                scale_color_manual(values = diet_colors) +
                scale_fill_manual(values = diet_colors) +
                scale_x_continuous(breaks = c(1, 7, 14, 21, 35, 56, 70)) +
                theme_bw(base_size = 10) +
                theme(legend.position = "none")
        }
    }
    
    if(length(traj_plots) > 0) {
        n_col <- min(3, length(traj_plots))
        traj_combined <- wrap_plots(traj_plots, ncol = n_col) +
            plot_annotation(
                title = "Trajectories of High-Confidence Differential Genera",
                caption = "Statistics from ANCOM-BC2 and MaAsLin3 shown"
            )
        
        ggsave("high_confidence_trajectories.pdf", traj_combined,
               width = 4 * n_col, height = 4 * ceiling(length(traj_plots)/n_col))
        print(traj_combined)
    }
}

# ============================================================================
# FINAL INTEGRATIVE REPORT
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("GENERATING FINAL INTEGRATIVE REPORT\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

sink("final_integrative_report.txt")

cat("========================================\n")
cat("FINAL INTEGRATIVE MICROBIOME ANALYSIS REPORT\n")
cat("========================================\n\n")
cat("Date:", date(), "\n")
cat("Methods: ANCOM-BC2 and MaAsLin3\n\n")

cat("EXECUTIVE SUMMARY\n")
cat("------------------\n")
cat("Two complementary differential abundance methods were used to\n")
cat("identify genera associated with HFD vs ND diet.\n\n")
cat("High-confidence genera (detected by both methods):", length(overlap_genera), "\n")
cat("Total differential genera identified:", length(all_sig_genera), "\n\n")

cat("HIGH-CONFIDENCE GENERA (BOTH METHODS)\n")
cat("--------------------------------------\n")
if(length(overlap_genera) > 0) {
    cat("The following genera were consistently identified by both methods:\n\n")
    for(g in overlap_genera) {
        ancom_idx <- which(ancom_results$taxon_clean == g)
        maaslin_idx <- which(maaslin_diet$taxon_clean == g)
        
        cat(paste0("  • ", g, "\n"))
        cat(paste0("    ANCOM-BC2: LFC = ", round(ancom_results$lfc[ancom_idx[1]], 3),
                   ", q = ", format(ancom_results$qval[ancom_idx[1]], scientific = TRUE, digits = 2), "\n"))
        cat(paste0("    MaAsLin3:  Coef = ", round(maaslin_diet$coef[maaslin_idx[1]], 3),
                   ", q = ", format(maaslin_diet$qval_individual[maaslin_idx[1]], scientific = TRUE, digits = 2), "\n\n"))
    }
}

cat("METHOD-SPECIFIC FINDINGS\n")
cat("------------------------\n")
cat("ANCOM-BC2 unique findings:", length(ancom_only), "\n")
if(length(ancom_only) > 0) {
    for(g in ancom_only) {
        cat(paste0("  • ", g, "\n"))
    }
}

cat("\nMaAsLin3 unique findings:", length(maaslin_only), "\n")
if(length(maaslin_only) > 0) {
    for(g in maaslin_only) {
        cat(paste0("  • ", g, "\n"))
    }
}

cat("\nRECOMMENDATIONS\n")
cat("----------------\n")
cat("1. High-confidence genera (both methods) should be prioritized for validation\n")
cat("2. Genera detected by only one method should be considered suggestive\n")
cat("3. Consider biological plausibility and literature support for final selection\n")

cat("\n========================================\n")
cat("All results saved in:", output_dir, "/\n")
cat("========================================\n")

sink()

# ============================================================================
# SAVE FINAL WORKSPACE
# ============================================================================

save.image("final_integrative_analysis.RData")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("FINAL INTEGRATIVE ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
cat("Generated files:\n")
cat("  - method_comparison_venn.pdf\n")
cat("  - method_effect_concordance.pdf\n")
cat("  - combined_significant_genera_heatmap.pdf\n")
cat("  - high_confidence_trajectories.pdf\n")
cat("  - high_confidence_genera.csv\n")
cat("  - method_effect_comparison.csv\n")
cat("  - final_integrative_report.txt\n")
cat("  - final_integrative_analysis.RData\n")

cat("\nSummary:\n")
cat("  Total differential genera:", length(all_sig_genera), "\n")
cat("  High-confidence (both methods):", length(overlap_genera), "\n")
cat("  ANCOM-BC2 only:", length(ancom_only), "\n")
cat("  MaAsLin3 only:", length(maaslin_only), "\n")