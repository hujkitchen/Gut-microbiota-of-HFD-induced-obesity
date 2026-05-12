#=============================================================================
# ANALYZE AND VISUALIZE MaAsLin3 RESULTS
#=============================================================================

#load required packages

library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(patchwork)


# ============================================================================
# LOAD AND PROCESS MaAsLin3 RESULTS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYZING MaAsLin3 RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Read MaAsLin3 results
maaslin_results_file <- "maaslin3_diet_effect/significant_results.tsv"

if (file.exists(maaslin_results_file)) {
    # Read the results
    sig_results_full <- read.delim(maaslin_results_file, stringsAsFactors = FALSE)
    
    cat("MaAsLin3 results loaded successfully!\n")
    cat("Total significant associations:", nrow(sig_results_full), "\n")
    cat("Columns:", paste(colnames(sig_results_full), collapse = ", "), "\n\n")
    
    # Check structure
    cat("Metadata variables tested:\n")
    print(table(sig_results_full$metadata))
    
    cat("\nFeatures with significant associations per metadata:\n")
    print(table(sig_results_full$metadata, sig_results_full$name))
    
    # Display first few rows
    cat("\nFirst 10 rows:\n")
    print(head(sig_results_full, 10))
    
    # ========================================================================
    # SEPARATE RESULTS BY METADATA VARIABLE
    # ========================================================================
    
    # Diet effect (group1)
    diet_results <- sig_results_full %>%
        filter(metadata == "group1") %>%
        arrange(qval_individual)
    
    cat("\n", paste(rep("-", 60), collapse = ""), "\n")
    cat("DIET EFFECT RESULTS (HFD vs ND)\n")
    cat(paste(rep("-", 60), collapse = ""), "\n\n")
    cat("Genera with significant diet effect:", nrow(diet_results), "\n")
    
    if (nrow(diet_results) > 0) {
        cat("\nResults summary:\n")
        print(diet_results[, c("feature", "coef", "pval_individual", "qval_individual", "model")])
        
        # Count by direction
        up_hfd <- sum(diet_results$coef > 0)
        down_hfd <- sum(diet_results$coef < 0)
        cat("\nGenera increased in HFD:", up_hfd, "\n")
        cat("Genera decreased in HFD:", down_hfd, "\n")
    }
    
    # Day effect
    day_results <- sig_results_full %>%
        filter(metadata == "day_scaled") %>%
        arrange(qval_individual)
    
    cat("\n", paste(rep("-", 60), collapse = ""), "\n")
    cat("DAY EFFECT RESULTS\n")
    cat(paste(rep("-", 60), collapse = ""), "\n\n")
    cat("Genera with significant temporal trend:", nrow(day_results), "\n")
    
    if (nrow(day_results) > 0) {
        cat("\nResults summary:\n")
        print(day_results[, c("feature", "coef", "pval_individual", "qval_individual")])
    }
    
    # Weight effect
    weight_results <- sig_results_full %>%
        filter(metadata == "weight_scaled") %>%
        arrange(qval_individual)
    
    cat("\n", paste(rep("-", 60), collapse = ""), "\n")
    cat("WEIGHT EFFECT RESULTS\n")
    cat(paste(rep("-", 60), collapse = ""), "\n\n")
    cat("Genera significantly associated with weight:", nrow(weight_results), "\n")
    
    if (nrow(weight_results) > 0) {
        cat("\nResults summary:\n")
        print(weight_results[, c("feature", "coef", "pval_individual", "qval_individual")])
    }
    
    # ========================================================================
    # VOLCANO PLOT FOR DIET EFFECT
    # ========================================================================
    
    if (nrow(diet_results) > 0) {
        cat("\n", paste(rep("=", 80), collapse = ""), "\n")
        cat("VOLCANO PLOT - DIET EFFECT\n")
        cat(paste(rep("=", 80), collapse = ""), "\n\n")
        
        # Prepare data
        volcano_df <- diet_results %>%
            mutate(
                significant = case_when(
                    qval_individual < 0.05 & coef > 0 ~ "Up in HFD",
                    qval_individual < 0.05 & coef < 0 ~ "Down in HFD",
                    TRUE ~ "NS"
                ),
                neg_log10_p = -log10(pval_individual),
                # Clean up feature names
                feature_clean = gsub("^g__", "", feature)
            )
        
        # Volcano plot
        p_volcano <- ggplot(volcano_df, aes(x = coef, y = neg_log10_p)) +
            geom_point(aes(color = significant), alpha = 0.7, size = 3) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 1) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
            geom_text_repel(
                data = subset(volcano_df, qval_individual < 0.01),
                aes(label = feature_clean),
                size = 4,
                max.overlaps = 20,
                box.padding = 0.5,
                fontface = "italic"
            ) +
            labs(
                x = "Coefficient (HFD vs ND)",
                y = "-log10(P-value)",
                title = "Differential Abundance: Diet Effect (MaAsLin3)",
                subtitle = paste0(
                    "HFD vs ND: ", 
                    sum(volcano_df$significant == "Up in HFD"), " increased, ",
                    sum(volcano_df$significant == "Down in HFD"), " decreased"
                )
            ) +
            scale_color_manual(
                values = c(
                    "Up in HFD" = "#A23B72", 
                    "Down in HFD" = "#2E86AB",
                    "NS" = "grey70"
                ),
                name = "Direction"
            ) +
            theme_bw(base_size = 12) +
            theme(
                legend.position = "bottom",
                plot.title = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(size = 11)
            )
        
        ggsave("maaslin3_volcano_diet.pdf", p_volcano, width = 10, height = 8)
        print(p_volcano)
        
        # Save volcano plot data
        write.csv(volcano_df, "volcano_plot_data.csv", row.names = FALSE)
    }
    
    # ========================================================================
    # HEATMAP OF ALL SIGNIFICANT GENERA
    # ========================================================================
    
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("HEATMAP OF SIGNIFICANT GENERA\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
    
    # Get unique significant genera
    sig_genera <- unique(sig_results_full$feature)
    
    if (length(sig_genera) > 0 && exists("genus_rel")) {
        # Extract relative abundances for significant genera
        sig_genera_present <- sig_genera[sig_genera %in% rownames(genus_rel)]
        
        if (length(sig_genera_present) > 0) {
            heatmap_data <- genus_rel[sig_genera_present, , drop = FALSE]
            
            # Log transform
            heatmap_data_log <- log10(heatmap_data + 1e-5)
            
            # Clean up feature names
            rownames(heatmap_data_log) <- gsub("^g__", "", rownames(heatmap_data_log))
            
            # Create annotation
            anno_col <- data.frame(
                Diet = metadata_clean$group1,
                Phase = metadata_clean$Phase,
                row.names = colnames(heatmap_data_log)
            )
            
            anno_colors <- list(
                Diet = diet_colors,
                Phase = phase_colors
            )
            
            # Heatmap
            pdf("maaslin3_significant_genera_heatmap.pdf", width = 14, height = max(8, length(sig_genera_present) * 0.3))
            
            if (nrow(heatmap_data_log) > 1) {
                pheatmap(heatmap_data_log,
                         scale = "row",
                         clustering_distance_rows = "correlation",
                         clustering_distance_cols = "euclidean",
                         annotation_col = anno_col,
                         annotation_colors = anno_colors,
                         show_rownames = TRUE,
                         fontsize_row = 9,
                         main = paste("Significant Genera from MaAsLin3 (n =", length(sig_genera_present), ")"),
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
            } else {
                cat("Only one significant genus, skipping heatmap\n")
            }
            
            dev.off()
            cat("Heatmap saved to maaslin3_significant_genera_heatmap.pdf\n")
        }
    }
    
    # ========================================================================
    # DOT PLOT FOR DIET EFFECT
    # ========================================================================
    
    if (nrow(diet_results) > 0) {
        cat("\n", paste(rep("=", 80), collapse = ""), "\n")
        cat("DOT PLOT - DIET EFFECT\n")
        cat(paste(rep("=", 80), collapse = ""), "\n\n")
        
        # Prepare data for dot plot
        dot_data <- diet_results %>%
            mutate(
                feature_clean = gsub("^g__", "", feature),
                neg_log10_q = -log10(qval_individual),
                direction = ifelse(coef > 0, "Up in HFD", "Down in HFD")
            ) %>%
            arrange(coef)
        
        # Dot plot
        p_dot <- ggplot(dot_data, aes(x = coef, y = reorder(feature_clean, coef))) +
            geom_point(aes(size = neg_log10_q, color = direction), alpha = 0.8) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
            labs(
                x = "Coefficient (HFD vs ND)",
                y = "",
                title = "Diet Effect on Genera Abundance",
                subtitle = "Point size indicates significance level"
            ) +
            scale_color_manual(values = c("Up in HFD" = "#A23B72", "Down in HFD" = "#2E86AB")) +
            scale_size_continuous(name = "-log10(Q-value)") +
            theme_bw(base_size = 12) +
            theme(
                axis.text.y = element_text(face = "italic"),
                legend.position = "bottom"
            )
        
        ggsave("maaslin3_dotplot_diet.pdf", p_dot, width = 10, height = max(6, nrow(diet_results) * 0.3))
        print(p_dot)
    }
    
    # ========================================================================
    # TRAJECTORIES OF TOP DIET-ASSOCIATED GENERA
    # ========================================================================
    
    if (nrow(diet_results) > 0 && exists("genus_rel")) {
        cat("\n", paste(rep("=", 80), collapse = ""), "\n")
        cat("TRAJECTORIES OF TOP GENERA\n")
        cat(paste(rep("=", 80), collapse = ""), "\n\n")
        
        # Get top 9 genera by significance
        top_n <- min(9, nrow(diet_results))
        top_genera <- diet_results %>%
            arrange(qval_individual) %>%
            slice(1:top_n) %>%
            pull(feature)
        
        # Prepare trajectory data
        traj_list <- list()
        for (g in top_genera) {
            if (g %in% rownames(genus_rel)) {
                genus_clean <- gsub("^g__", "", g)
                
                traj_df <- data.frame(
                    Abundance = genus_rel[g, ],
                    Day = metadata_clean$Day,
                    Diet = metadata_clean$group1,
                    Phase = metadata_clean$Phase
                )
                
                # Summary statistics
                traj_summary <- traj_df %>%
                    group_by(Diet, Day) %>%
                    summarise(
                        Mean = mean(Abundance),
                        SE = sd(Abundance) / sqrt(n()),
                        .groups = 'drop'
                    )
                
                # Get q-value for annotation
                qval_g <- diet_results$qval_individual[diet_results$feature == g]
                
                traj_list[[genus_clean]] <- ggplot(traj_summary, aes(x = Day, y = Mean, color = Diet)) +
                    geom_line(size = 1.2) +
                    geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE, fill = Diet), 
                                alpha = 0.2, color = NA) +
                    geom_point(data = traj_df, aes(x = Day, y = Abundance, color = Diet),
                               alpha = 0.2, size = 0.5) +
                    labs(
                        x = "Day", 
                        y = "Relative Abundance",
                        title = bquote(italic(.(genus_clean)) ~ "(q =" ~ .(format(qval_g, scientific = TRUE, digits = 2)) ~ ")")
                    ) +
                    scale_color_manual(values = diet_colors) +
                    scale_fill_manual(values = diet_colors) +
                    scale_x_continuous(breaks = c(1, 7, 14, 21, 35, 56, 70)) +
                    theme_bw(base_size = 10) +
                    theme(legend.position = "none")
            }
        }
        
        if (length(traj_list) > 0) {
            traj_combined <- wrap_plots(traj_list, ncol = 3) +
                plot_annotation(
                    title = "Trajectories of Top Diet-Associated Genera",
                    theme = theme(plot.title = element_text(face = "bold", size = 14))
                )
            
            ggsave("maaslin3_top_genera_trajectories.pdf", traj_combined, 
                   width = 14, height = ceiling(length(traj_list)/3) * 4)
            print(traj_combined)
        }
    }
    
    # ========================================================================
    # DAY EFFECT VISUALIZATION
    # ========================================================================
    
    if (nrow(day_results) > 0) {
        cat("\n", paste(rep("=", 80), collapse = ""), "\n")
        cat("TEMPORAL TRENDS VISUALIZATION\n")
        cat(paste(rep("=", 80), collapse = ""), "\n\n")
        
        # Get top genera with temporal trends
        top_day_genera <- day_results %>%
            arrange(qval_individual) %>%
            slice(1:min(6, n())) %>%
            mutate(feature_clean = gsub("^g__", "", feature))
        
        # Create scatter plots with trend lines
        day_plot_list <- list()
        for (i in 1:nrow(top_day_genera)) {
            g <- top_day_genera$feature[i]
            g_clean <- top_day_genera$feature_clean[i]
            
            if (g %in% rownames(genus_rel)) {
                day_df <- data.frame(
                    Abundance = genus_rel[g, ],
                    Day = metadata_clean$Day,
                    Diet = metadata_clean$group1
                )
                
                day_plot_list[[g_clean]] <- ggplot(day_df, aes(x = Day, y = Abundance, color = Diet)) +
                    geom_point(alpha = 0.5, size = 1.5) +
                    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
                    facet_wrap(~ Diet) +
                    labs(
                        x = "Day", 
                        y = "Relative Abundance",
                        title = bquote(italic(.(g_clean)))
                    ) +
                    scale_color_manual(values = diet_colors) +
                    theme_bw(base_size = 10) +
                    theme(legend.position = "none")
            }
        }
        
        if (length(day_plot_list) > 0) {
            day_combined <- wrap_plots(day_plot_list, ncol = 2)
            ggsave("maaslin3_temporal_trends.pdf", day_combined, width = 12, height = ceiling(length(day_plot_list)/2) * 4)
            print(day_combined)
        }
    }
    
    # ========================================================================
    # WEIGHT ASSOCIATION VISUALIZATION
    # ========================================================================
    
    if (nrow(weight_results) > 0) {
        cat("\n", paste(rep("=", 80), collapse = ""), "\n")
        cat("WEIGHT ASSOCIATION VISUALIZATION\n")
        cat(paste(rep("=", 80), collapse = ""), "\n\n")
        
        # Get top weight-associated genera
        top_weight_genera <- weight_results %>%
            arrange(qval_individual) %>%
            slice(1:min(6, n())) %>%
            mutate(feature_clean = gsub("^g__", "", feature))
        
        weight_plot_list <- list()
        for (i in 1:nrow(top_weight_genera)) {
            g <- top_weight_genera$feature[i]
            g_clean <- top_weight_genera$feature_clean[i]
            
            if (g %in% rownames(genus_rel)) {
                weight_df <- data.frame(
                    Abundance = genus_rel[g, ],
                    Weight = metadata_clean$weight,
                    Day = metadata_clean$Day,
                    Diet = metadata_clean$group1
                )
                
                weight_plot_list[[g_clean]] <- ggplot(weight_df, aes(x = Abundance, y = Weight, color = Diet)) +
                    geom_point(alpha = 0.6, size = 2) +
                    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
                    labs(
                        x = paste0("Relative Abundance of ", g_clean),
                        y = "Weight (g)",
                        title = bquote(italic(.(g_clean)))
                    ) +
                    scale_color_manual(values = diet_colors) +
                    theme_bw(base_size = 10) +
                    theme(legend.position = "bottom")
            }
        }
        
        if (length(weight_plot_list) > 0) {
            weight_combined <- wrap_plots(weight_plot_list, ncol = 2)
            ggsave("maaslin3_weight_associations.pdf", weight_combined, 
                   width = 12, height = ceiling(length(weight_plot_list)/2) * 4)
            print(weight_combined)
        }
    }
    
    # ========================================================================
    # SUMMARY TABLE
    # ========================================================================
    
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("MAASLIN3 RESULTS SUMMARY\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
    
    summary_table <- sig_results_full %>%
        group_by(metadata, name) %>%
        summarise(
            n_genera = n(),
            mean_coef = mean(coef),
            min_qval = min(qval_individual),
            .groups = 'drop'
        ) %>%
        arrange(min_qval)
    
    print(summary_table)
    write.csv(summary_table, "maaslin3_results_summary.csv", row.names = FALSE)
    
    # Export all significant results
    write.csv(sig_results_full, "maaslin3_all_significant_results.csv", row.names = FALSE)
    
    # ========================================================================
    # FINAL REPORT
    # ========================================================================
    
    sink("maaslin3_analysis_report.txt")
    
    cat("========================================\n")
    cat("MaAsLin3 DIFFERENTIAL ABUNDANCE ANALYSIS REPORT\n")
    cat("========================================\n\n")
    cat("Date:", date(), "\n\n")
    
    cat("1. OVERVIEW\n")
    cat("-----------\n")
    cat("Total significant associations:", nrow(sig_results_full), "\n")
    cat("Unique genera with significant associations:", length(unique(sig_results_full$feature)), "\n\n")
    
    cat("2. DIET EFFECT (HFD vs ND)\n")
    cat("--------------------------\n")
    if (nrow(diet_results) > 0) {
        cat("Genera significantly affected by diet:", nrow(diet_results), "\n")
        cat("  - Increased in HFD:", sum(diet_results$coef > 0), "\n")
        cat("  - Decreased in HFD:", sum(diet_results$coef < 0), "\n\n")
        cat("Top 10 genera:\n")
        diet_top10 <- head(diet_results[, c("feature", "coef", "qval_individual", "model")], 10)
        print(diet_top10)
    } else {
        cat("No genera significantly affected by diet.\n")
    }
    
    cat("\n3. TEMPORAL TRENDS\n")
    cat("-----------------\n")
    if (nrow(day_results) > 0) {
        cat("Genera with significant temporal trends:", nrow(day_results), "\n")
        cat("  - Positive trend:", sum(day_results$coef > 0), "\n")
        cat("  - Negative trend:", sum(day_results$coef < 0), "\n\n")
        cat("Top 10 genera:\n")
        day_top10 <- head(day_results[, c("feature", "coef", "qval_individual")], 10)
        print(day_top10)
    } else {
        cat("No genera with significant temporal trends.\n")
    }
    
    cat("\n4. WEIGHT ASSOCIATIONS\n")
    cat("---------------------\n")
    if (nrow(weight_results) > 0) {
        cat("Genera significantly associated with weight:", nrow(weight_results), "\n")
        cat("  - Positive association:", sum(weight_results$coef > 0), "\n")
        cat("  - Negative association:", sum(weight_results$coef < 0), "\n\n")
        cat("Top 10 genera:\n")
        weight_top10 <- head(weight_results[, c("feature", "coef", "qval_individual")], 10)
        print(weight_top10)
    } else {
        cat("No genera significantly associated with weight.\n")
    }
    
    cat("\n========================================\n")
    cat("All results saved in:", output_dir, "/\n")
    cat("========================================\n")
    
    sink()
    
    cat("\nMaAsLin3 analysis complete!\n")
    cat("Report saved to: maaslin3_analysis_report.txt\n")
    
} else {
    cat("ERROR: MaAsLin3 results file not found at:", maaslin_results_file, "\n")
    cat("Please check the file path.\n")
}

# ============================================================================
# SAVE FINAL WORKSPACE
# ============================================================================

save.image("complete_analysis_with_maaslin3.RData")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ALL ANALYSES COMPLETE!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
cat("Generated files:\n")
cat("  - maaslin3_volcano_diet.pdf\n")
cat("  - maaslin3_dotplot_diet.pdf\n")
cat("  - maaslin3_top_genera_trajectories.pdf\n")
cat("  - maaslin3_temporal_trends.pdf\n")
cat("  - maaslin3_weight_associations.pdf\n")
cat("  - maaslin3_significant_genera_heatmap.pdf\n")
cat("  - maaslin3_all_significant_results.csv\n")
cat("  - maaslin3_results_summary.csv\n")
cat("  - maaslin3_analysis_report.txt\n")
cat("  - complete_analysis_with_maaslin3.RData\n")