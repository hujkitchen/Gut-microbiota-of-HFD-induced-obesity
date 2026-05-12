#=============================================================================
# COMPREHENSIVE PLOTS FOR THE FIVE HIGH-CONFIDENCE GENERA
#=============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpubr)
library(viridis)

# Define the five genera
five_genera <- c("Lactococcus", "Duncaniella", "Taurinivorans", "Lepagella", "CAG-873")

# Colors
diet_colors <- c("ND" = "#2E86AB", "HFD" = "#A23B72")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("PLOTTING FIVE HIGH-CONFIDENCE DIFFERENTIAL GENERA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ============================================================================
# CHECK WHICH GENERA ARE IN THE DATA
# ============================================================================

cat("Checking availability of genera in data...\n")
for(g in five_genera) {
    if(g %in% rownames(genus_rel)) {
        cat("  ✓", g, "- Found in genus_rel\n")
    } else {
        # Try alternative names
        g_alt <- paste0("g__", g)
        if(g_alt %in% rownames(genus_rel)) {
            cat("  ✓", g, "- Found as", g_alt, "\n")
        } else {
            # Try case-insensitive
            match_idx <- grep(g, rownames(genus_rel), ignore.case = TRUE)
            if(length(match_idx) > 0) {
                cat("  ✓", g, "- Found as", rownames(genus_rel)[match_idx[1]], "\n")
            } else {
                cat("  ✗", g, "- NOT FOUND\n")
            }
        }
    }
}

# ============================================================================
# EXTRACT DATA FOR FIVE GENERA
# ============================================================================

# Create a function to get genus data
get_genus_data <- function(genus_name) {
    # Try exact match
    if(genus_name %in% rownames(genus_rel)) {
        return(genus_rel[genus_name, ])
    }
    # Try with g__ prefix
    g_alt <- paste0("g__", genus_name)
    if(g_alt %in% rownames(genus_rel)) {
        return(genus_rel[g_alt, ])
    }
    # Try case-insensitive
    match_idx <- grep(genus_name, rownames(genus_rel), ignore.case = TRUE)
    if(length(match_idx) > 0) {
        return(genus_rel[match_idx[1], ])
    }
    return(NULL)
}

# Build combined dataframe
five_genera_data <- data.frame()
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
        five_genera_data <- rbind(five_genera_data, temp_df)
    }
}

# Add log-transformed abundance
five_genera_data$LogAbundance <- log10(five_genera_data$Abundance + 1e-6)

cat("\nData extracted for", length(unique(five_genera_data$Genus)), "genera\n")
cat("Total observations:", nrow(five_genera_data), "\n")

# ============================================================================
# CALCULATE STATISTICS FOR EACH GENUS
# ============================================================================

genus_stats <- five_genera_data %>%
    group_by(Genus, Diet) %>%
    summarise(
        Mean = mean(Abundance),
        SD = sd(Abundance),
        SEM = sd(Abundance) / sqrt(n()),
        Median = median(Abundance),
        .groups = 'drop'
    )

# Calculate fold changes
genus_fc <- genus_stats %>%
    select(Genus, Diet, Mean) %>%
    pivot_wider(names_from = Diet, values_from = Mean) %>%
    mutate(
        FC_HFD_ND = HFD / ND,
        Log2FC = log2(HFD / ND),
        Direction = ifelse(Log2FC > 0, "↑ Up in HFD", "↓ Down in HFD")
    )

cat("\nGenus statistics:\n")
print(genus_fc)

# Wilcoxon test per genus
wilcox_results <- data.frame()
for(g in five_genera) {
    g_data <- five_genera_data %>% filter(Genus == g)
    if(nrow(g_data) > 0 && length(unique(g_data$Diet)) == 2) {
        nd_vals <- g_data$Abundance[g_data$Diet == "ND"]
        hfd_vals <- g_data$Abundance[g_data$Diet == "HFD"]
        test <- wilcox.test(hfd_vals, nd_vals)
        wilcox_results <- rbind(wilcox_results, data.frame(
            Genus = g,
            W = test$statistic,
            p_value = test$p.value,
            stringsAsFactors = FALSE
        ))
    }
}
wilcox_results$p_adj <- p.adjust(wilcox_results$p_value, method = "BH")
wilcox_results$significance <- ifelse(wilcox_results$p_adj < 0.001, "***",
                                ifelse(wilcox_results$p_adj < 0.01, "**",
                                ifelse(wilcox_results$p_adj < 0.05, "*", "ns")))

cat("\nWilcoxon test results:\n")
print(wilcox_results)

# ============================================================================
# CALCULATE SUMMARY FOR TRAJECTORIES
# ============================================================================

traj_summary <- five_genera_data %>%
    group_by(Genus, Diet, Day) %>%
    summarise(
        Mean = mean(Abundance),
        SE = sd(Abundance) / sqrt(n()),
        n = n(),
        .groups = 'drop'
    )

# Smooth trend using loess
traj_smooth <- five_genera_data %>%
    group_by(Genus, Diet) %>%
    do({
        model <- loess(Abundance ~ Day, data = ., span = 0.5)
        day_seq <- seq(min(.$Day), max(.$Day), length.out = 100)
        pred <- predict(model, newdata = data.frame(Day = day_seq), se = TRUE)
        data.frame(
            Day = day_seq,
            Smooth = pred$fit,
            SE_Smooth = pred$se.fit,
            Genus = unique(.$Genus),
            Diet = unique(.$Diet)
        )
    })

# ============================================================================
# PLOT 1: BOXPLOT WITH JITTER - DIET COMPARISON
# ============================================================================

p_boxplot <- ggplot(five_genera_data, aes(x = Diet, y = Abundance, fill = Diet)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = Diet), width = 0.15, alpha = 0.4, size = 1) +
    facet_wrap(~ Genus, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = diet_colors, guide = "none") +
    scale_color_manual(values = diet_colors, guide = "none") +
    scale_y_log10(labels = scales::scientific) +
    labs(
        title = "Abundance Distribution by Diet",
        subtitle = "Five high-confidence differential genera",
        x = "",
        y = "Relative Abundance (log10 scale)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        strip.text = element_text(face = "italic", size = 11),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(face = "bold")
    ) +
    # Add significance from Wilcoxon test
    geom_text(
        data = wilcox_results,
        aes(x = 1.5, y = Inf, label = significance),
        inherit.aes = FALSE,
        vjust = 1.5,
        size = 6,
        color = "black"
    )

ggsave("five_genera_boxplot.pdf", p_boxplot, width = 14, height = 5)
print(p_boxplot)

# ============================================================================
# PLOT 2: VIOLIN PLOT WITH BOXPLOT
# ============================================================================

p_violin <- ggplot(five_genera_data, aes(x = Diet, y = Abundance, fill = Diet)) +
    geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75), trim = TRUE) +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
    facet_wrap(~ Genus, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = diet_colors, guide = "none") +
    scale_y_log10(labels = scales::scientific) +
    labs(
        title = "Abundance Distribution: Violin Plot",
        subtitle = "With quartiles and median",
        x = "",
        y = "Relative Abundance (log10 scale)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        strip.text = element_text(face = "italic", size = 11),
        plot.title = element_text(face = "bold")
    )

ggsave("five_genera_violin.pdf", p_violin, width = 14, height = 5)
print(p_violin)

# ============================================================================
# PLOT 3: TIME SERIES TRAJECTORIES WITH LOESS SMOOTHING
# ============================================================================

p_trajectory <- ggplot() +
    # Individual points
    geom_point(data = five_genera_data,
               aes(x = Day, y = Abundance, color = Diet),
               alpha = 0.2, size = 0.8) +
    # Summary line with error ribbon
    geom_ribbon(data = traj_summary,
                aes(x = Day, ymin = Mean - SE, ymax = Mean + SE, fill = Diet),
                alpha = 0.2) +
    geom_line(data = traj_summary,
              aes(x = Day, y = Mean, color = Diet),
              linewidth = 1.2) +
    # LOESS smooth
    geom_line(data = traj_smooth,
              aes(x = Day, y = Smooth, color = Diet),
              linewidth = 0.8, linetype = "dashed", alpha = 0.7) +
    facet_wrap(~ Genus, scales = "free_y", ncol = 3) +
    scale_color_manual(values = diet_colors, name = "Diet") +
    scale_fill_manual(values = diet_colors, name = "Diet") +
    scale_x_continuous(breaks = c(1, 7, 14, 21, 35, 56, 70)) +
    labs(
        title = "Temporal Trajectories of Differential Genera",
        subtitle = "Solid line: daily mean ± SE | Dashed line: LOESS smooth (span = 0.5)",
        x = "Day",
        y = "Relative Abundance"
    ) +
    theme_bw(base_size = 12) +
    theme(
        strip.text = element_text(face = "italic", size = 11),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
    )

ggsave("five_genera_trajectories.pdf", p_trajectory, width = 14, height = 10)
print(p_trajectory)

# ============================================================================
# PLOT 4: PHASE BOXPLOT
# ============================================================================

p_phase <- ggplot(five_genera_data, aes(x = Phase, y = Abundance, fill = Diet)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8), outlier.shape = NA) +
    geom_jitter(aes(color = Diet), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                alpha = 0.2, size = 0.5) +
    facet_wrap(~ Genus, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = diet_colors, name = "Diet") +
    scale_color_manual(values = diet_colors, guide = "none") +
    scale_y_log10(labels = scales::scientific) +
    labs(
        title = "Abundance by Study Phase",
        subtitle = "Acute (D1-7) | Early (D14-21) | Mid (D28-42) | Late (D56-70)",
        x = "",
        y = "Relative Abundance (log10 scale)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        strip.text = element_text(face = "italic", size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
    )

ggsave("five_genera_by_phase.pdf", p_phase, width = 14, height = 10)
print(p_phase)

# ============================================================================
# PLOT 5: HEATMAP-STYLE DOT PLOT
# ============================================================================

# Prepare data for heatmap
heatmap_data <- five_genera_data %>%
    group_by(Genus, Diet, Day) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = 'drop') %>%
    mutate(
        LogMean = log10(MeanAbundance + 1e-6),
        Day = factor(Day, levels = sort(unique(Day)))
    )

p_heatmap <- ggplot(heatmap_data, aes(x = Day, y = Genus, fill = LogMean)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_wrap(~ Diet, ncol = 1) +
    scale_fill_viridis(option = "plasma", name = "log10(Abundance)") +
    labs(
        title = "Temporal Abundance Heatmap",
        subtitle = "Five high-confidence differential genera",
        x = "Day",
        y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.y = element_text(face = "italic", size = 11),
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid = element_blank()
    )

ggsave("five_genera_heatmap.pdf", p_heatmap, width = 12, height = 6)
print(p_heatmap)

# ============================================================================
# PLOT 6: COMBINED EFFECT SIZE PLOT
# ============================================================================

# Get effect sizes from validation_results
effect_data <- validation_results %>%
    filter(genus %in% five_genera) %>%
    select(genus, true_log2FC, ancom_lfc_corrected, maaslin_coef) %>%
    pivot_longer(cols = c(true_log2FC, ancom_lfc_corrected, maaslin_coef),
                 names_to = "Method",
                 values_to = "EffectSize") %>%
    mutate(
        Method = case_when(
            Method == "true_log2FC" ~ "Actual Data",
            Method == "ancom_lfc_corrected" ~ "ANCOM-BC2 (Corrected)",
            Method == "maaslin_coef" ~ "MaAsLin3"
        ),
        Method = factor(Method, levels = c("Actual Data", "ANCOM-BC2 (Corrected)", "MaAsLin3"))
    )

p_effect <- ggplot(effect_data, aes(x = reorder(genus, EffectSize), y = EffectSize, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), alpha = 0.8, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(
        values = c("Actual Data" = "#4CAF50", 
                   "ANCOM-BC2 (Corrected)" = "#2E86AB", 
                   "MaAsLin3" = "#A23B72"),
        name = ""
    ) +
    labs(
        title = "Effect Size Comparison: HFD vs ND",
        subtitle = "Positive = higher in HFD | Negative = lower in HFD",
        x = "",
        y = "Effect Size"
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.y = element_text(face = "italic", size = 12),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
    ) +
    # Add value labels
    geom_text(
        aes(label = round(EffectSize, 2)),
        position = position_dodge(0.8),
        hjust = ifelse(effect_data$EffectSize > 0, -0.3, 1.3),
        size = 3.5
    )

ggsave("five_genera_effect_sizes.pdf", p_effect, width = 10, height = 6)
print(p_effect)

# ============================================================================
# PLOT 7: INDIVIDUAL DETAILED PLOTS FOR EACH GENUS
# ============================================================================

individual_plots <- list()

for(g in five_genera) {
    g_data <- five_genera_data %>% filter(Genus == g)
    
    if(nrow(g_data) > 0) {
        # Get stats
        fc_info <- genus_fc %>% filter(Genus == g)
        w_test <- wilcox_results %>% filter(Genus == g)
        
        # Create multi-panel plot for each genus
        p1 <- ggplot(g_data, aes(x = Diet, y = Abundance, fill = Diet)) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            geom_jitter(aes(color = Diet), width = 0.15, alpha = 0.5, size = 1.5) +
            scale_fill_manual(values = diet_colors) +
            scale_color_manual(values = diet_colors) +
            scale_y_log10() +
            labs(title = paste0("Abundance: ", g),
                 subtitle = paste0("Log2FC = ", round(fc_info$Log2FC, 2),
                                   " | p.adj = ", format(w_test$p_adj, scientific = TRUE, digits = 2),
                                   " ", w_test$significance),
                 x = "", y = "Relative Abundance") +
            theme_bw() +
            theme(legend.position = "none")
        
        p2 <- ggplot(g_data, aes(x = Day, y = Abundance, color = Diet)) +
            geom_point(alpha = 0.4, size = 1.5) +
            geom_smooth(aes(fill = Diet), method = "loess", se = TRUE, alpha = 0.2, span = 0.5) +
            scale_color_manual(values = diet_colors) +
            scale_fill_manual(values = diet_colors) +
            scale_x_continuous(breaks = c(1, 7, 14, 21, 35, 56, 70)) +
            labs(x = "Day", y = "Relative Abundance",
                 title = "Temporal Trend") +
            theme_bw() +
            theme(legend.position = "none")
        
        p3 <- ggplot(g_data, aes(x = Weight, y = Abundance, color = Diet)) +
            geom_point(alpha = 0.6, size = 2) +
            geom_smooth(aes(fill = Diet), method = "lm", se = TRUE, alpha = 0.2) +
            scale_color_manual(values = diet_colors) +
            scale_fill_manual(values = diet_colors) +
            labs(x = "Weight (g)", y = "Relative Abundance",
                 title = "Weight Association") +
            theme_bw() +
            theme(legend.position = "bottom")
        
        combined <- (p1 | p2 | p3) +
            plot_annotation(
                title = bquote(italic(.(g))),
                theme = theme(plot.title = element_text(size = 14, face = "bold.italic"))
            )
        
        individual_plots[[g]] <- combined
        
        ggsave(paste0("genus_", g, "_detailed.pdf"), combined, width = 14, height = 5)
    }
}

# Display combined plot
if(length(individual_plots) > 0) {
    all_individual <- wrap_plots(individual_plots, ncol = 1)
    ggsave("five_genera_detailed_combined.pdf", all_individual, 
           width = 16, height = 5 * length(individual_plots))
}

# ============================================================================
# PLOT 8: ABUNDANCE RANK PLOT
# ============================================================================

# Calculate mean abundance per genus
rank_data <- five_genera_data %>%
    group_by(Genus, Diet) %>%
    summarise(MeanAbund = mean(Abundance), .groups = 'drop') %>%
    group_by(Diet) %>%
    arrange(desc(MeanAbund)) %>%
    mutate(Rank = row_number())

p_rank <- ggplot(rank_data, aes(x = reorder(Genus, MeanAbund), y = MeanAbund, fill = Diet)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), alpha = 0.8, width = 0.6) +
    scale_fill_manual(values = diet_colors) +
    scale_y_log10(labels = scales::scientific) +
    coord_flip() +
    labs(
        title = "Mean Abundance Ranking",
        subtitle = "Comparison between ND and HFD",
        x = "",
        y = "Mean Relative Abundance (log10 scale)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.y = element_text(face = "italic", size = 12),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
    ) +
    geom_text(aes(label = format(MeanAbund, scientific = TRUE, digits = 2)),
              position = position_dodge(0.7),
              hjust = -0.2, size = 3)

ggsave("five_genera_abundance_rank.pdf", p_rank, width = 8, height = 5)
print(p_rank)

# ============================================================================
# COMBINED SUMMARY FIGURE
# ============================================================================

# Create a publication-ready combined figure
p_top <- p_boxplot + 
    labs(tag = "A") +
    theme(plot.tag = element_text(face = "bold", size = 14))

p_middle <- p_trajectory + 
    labs(tag = "B") +
    theme(plot.tag = element_text(face = "bold", size = 14))

p_bottom <- p_effect + 
    labs(tag = "C") +
    theme(plot.tag = element_text(face = "bold", size = 14))

combined_figure <- (p_top / p_middle / p_bottom) +
    plot_layout(heights = c(1, 1.5, 1))

ggsave("five_genera_summary_figure.pdf", combined_figure, 
       width = 16, height = 20)
print(combined_figure)

# ============================================================================
# SAVE DATA AND STATISTICS
# ============================================================================

# Save the processed data
write.csv(five_genera_data, "five_genera_abundance_data.csv", row.names = FALSE)
write.csv(genus_fc, "five_genera_fold_changes.csv", row.names = FALSE)
write.csv(wilcox_results, "five_genera_wilcoxon_test.csv", row.names = FALSE)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("PLOTS FOR FIVE GENERA COMPLETE!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
cat("Generated files:\n")
cat("  - five_genera_boxplot.pdf\n")
cat("  - five_genera_violin.pdf\n")
cat("  - five_genera_trajectories.pdf\n")
cat("  - five_genera_by_phase.pdf\n")
cat("  - five_genera_heatmap.pdf\n")
cat("  - five_genera_effect_sizes.pdf\n")
cat("  - five_genera_abundance_rank.pdf\n")
cat("  - five_genera_summary_figure.pdf\n")
cat("  - genus_*_detailed.pdf (individual detailed plots)\n")
cat("  - five_genera_detailed_combined.pdf\n")
cat("  - five_genera_abundance_data.csv\n")
cat("  - five_genera_fold_changes.csv\n")
cat("  - five_genera_wilcoxon_test.csv\n")

# Print summary of findings
cat("\nSummary of Five Genera:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
for(i in 1:nrow(genus_fc)) {
    w <- wilcox_results %>% filter(Genus == genus_fc$Genus[i])
    cat(sprintf("%-15s: Log2FC = %+6.2f | %s | p.adj = %s %s\n",
                genus_fc$Genus[i],
                genus_fc$Log2FC[i],
                genus_fc$Direction[i],
                ifelse(nrow(w) > 0, format(w$p_adj, scientific = TRUE, digits = 2), "NA"),
                ifelse(nrow(w) > 0, w$significance, "")))
}