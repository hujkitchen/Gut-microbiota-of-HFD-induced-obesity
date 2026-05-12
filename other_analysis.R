#=============================================================================
#LAG ANALYSIS AND EFFECT SIZE ANALYSIS
#=============================================================================


# The problem: five_genera has g__ prefix, data has clean names
# Let's fix five_genera once and for all

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("FIXING NAME MISMATCH\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("BEFORE:\n")
cat("  five_genera =", paste(five_genera, collapse=", "), "\n")
cat("  Data genera =", paste(unique(weight_cor_data$Genus), collapse=", "), "\n\n")

# Strip g__ prefix
five_genera <- gsub("^g__", "", five_genera)

cat("AFTER:\n")
cat("  five_genera =", paste(five_genera, collapse=", "), "\n")
cat("  Data genera =", paste(unique(weight_cor_data$Genus), collapse=", "), "\n\n")

# Verify match
matches <- five_genera %in% weight_cor_data$Genus
cat("Match status:\n")
for(i in 1:length(five_genera)) {
    cat(sprintf("  %s: %s\n", five_genera[i], ifelse(matches[i], "✓ FOUND", "✗ NOT FOUND")))
}

# ============================================================================
# NOW RUN LAG ANALYSIS WITH FIXED NAMES
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("LAG ANALYSIS - CROSS-CORRELATION \n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

lag_all_results <- data.frame()

for(g in five_genera) {
    cat("Processing:", g, "... ")
    
    g_data <- weight_cor_data %>% 
        filter(Genus == g) %>%
        arrange(Cage, Day)
    
    n_cages <- length(unique(g_data$Cage))
    cat(n_cages, "cages ... ")
    
    genus_lag <- data.frame()
    
    for(cage_id in unique(g_data$Cage)) {
        cage_data <- g_data %>% 
            filter(Cage == cage_id) %>%
            arrange(Day)
        
        if(nrow(cage_data) >= 6) {
            tryCatch({
                cc <- ccf(cage_data$Abundance, cage_data$Weight, 
                          lag.max = 4, plot = FALSE)
                
                genus_lag <- rbind(genus_lag, data.frame(
                    Lag = as.numeric(cc$lag),
                    Correlation = as.numeric(cc$acf),
                    Cage = cage_id,
                    Genus = g,
                    Diet = unique(cage_data$Diet),
                    stringsAsFactors = FALSE
                ))
            }, error = function(e) {})
        }
    }
    
    if(nrow(genus_lag) > 0) {
        lag_all_results <- rbind(lag_all_results, genus_lag)
        cat("OK (", nrow(genus_lag), " rows)\n", sep="")
    } else {
        cat("No results\n")
    }
}

# ============================================================================
# SUMMARIZE AND PLOT LAG RESULTS
# ============================================================================

if(nrow(lag_all_results) > 0) {
    
    lag_summary <- lag_all_results %>%
        group_by(Genus, Lag) %>%
        summarise(
            MeanCorrelation = mean(Correlation, na.rm = TRUE),
            SE = sd(Correlation, na.rm = TRUE) / sqrt(n()),
            N = n(),
            .groups = 'drop'
        )
    
    # Best lag per genus
    cat("\nOptimal lag (maximum absolute correlation):\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    for(g in unique(lag_summary$Genus)) {
        g_lag <- lag_summary %>% 
            filter(Genus == g, Lag != 0) %>%
            slice_max(order_by = abs(MeanCorrelation), n = 1)
        
        if(nrow(g_lag) > 0) {
            direction <- ifelse(g_lag$MeanCorrelation > 0, "positive", "negative")
            interpretation <- ifelse(g_lag$Lag < 0, 
                                     "Abundance PRECEDES weight change",
                                     "Weight change PRECEDES abundance shift")
            cat(sprintf("  %-20s: Lag = %+d | r = %+.3f (%s) | %s\n",
                        g, g_lag$Lag, g_lag$MeanCorrelation, direction, interpretation))
        }
    }
    
    # Plot
    genus_colors <- c(
        "Lactococcus" = "#E41A1C",
        "Duncaniella" = "#377EB8", 
        "Taurinivorans" = "#4DAF4A",
        "Lepagella" = "#984EA3",
        "CAG-873" = "#FF7F00"
    )
    
    p_lag <- ggplot(lag_summary, aes(x = Lag, y = MeanCorrelation, 
                                      color = Genus, fill = Genus)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.5) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        geom_ribbon(aes(ymin = MeanCorrelation - SE, ymax = MeanCorrelation + SE),
                    alpha = 0.12, color = NA) +
        facet_wrap(~ Genus, ncol = 3) +
        scale_color_manual(values = genus_colors, guide = "none") +
        scale_fill_manual(values = genus_colors, guide = "none") +
        scale_x_continuous(breaks = seq(-4, 4, by = 1)) +
        labs(
            x = "Lag (timepoints)",
            y = "Cross-Correlation (Abundance vs Weight)",
            title = "Temporal Lag Analysis",
            subtitle = "Negative lag → Abundance changes BEFORE weight | Positive lag → Weight changes BEFORE abundance"
        ) +
        theme_bw(base_size = 12) +
        theme(
            strip.text = element_text(face = "italic", size = 11),
            plot.title = element_text(face = "bold")
        )
    
    ggsave("lag_cross_correlation.pdf", p_lag, width = 13, height = 9)
    print(p_lag)
    
    write.csv(lag_summary, "lag_analysis_results.csv", row.names = FALSE)
    save(lag_all_results, lag_summary, file = "lag_analysis_results.RData")
    
    cat("\nLag analysis complete!\n")
    cat("  - lag_cross_correlation.pdf\n")
    cat("  - lag_analysis_results.csv\n")
    
} else {
    cat("\nERROR: Still no lag results. Check if abundance varies within cages.\n")
    for(g in five_genera) {
        g_data <- weight_cor_data %>% filter(Genus == g)
        cat(sprintf("  %s: non-zero values = %d/%d\n", 
                    g, sum(g_data$Abundance > 0), nrow(g_data)))
    }
}

# ============================================================================
# ALSO RUN EFFECT SIZE ANALYSIS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("EFFECT SIZE ANALYSIS (FIXED NAMES)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

if (!require("effsize", quietly = TRUE)) {
    install.packages("effsize")
    library(effsize)
}

effect_sizes <- data.frame()

for(g in five_genera) {
    g_data <- weight_cor_data %>% filter(Genus == g)
    
    nd_vals <- g_data$Abundance[g_data$Diet == "ND"]
    hfd_vals <- g_data$Abundance[g_data$Diet == "HFD"]
    nd_vals <- nd_vals[!is.na(nd_vals)]
    hfd_vals <- hfd_vals[!is.na(hfd_vals)]
    
    if(length(nd_vals) >= 3 && length(hfd_vals) >= 3) {
        cd <- cohen.d(hfd_vals, nd_vals, pooled = TRUE, paired = FALSE, na.rm = TRUE)
        
        effect_sizes <- rbind(effect_sizes, data.frame(
            Genus = g,
            Cohens_d = round(cd$estimate, 4),
            Magnitude = cd$magnitude,
            ND_mean = mean(nd_vals), HFD_mean = mean(hfd_vals),
            ND_sd = sd(nd_vals), HFD_sd = sd(hfd_vals),
            N_ND = length(nd_vals), N_HFD = length(hfd_vals),
            Log2FC = log2((mean(hfd_vals) + 1e-10) / (mean(nd_vals) + 1e-10)),
            stringsAsFactors = FALSE
        ))
    }
}

if(nrow(effect_sizes) > 0) {
    print(effect_sizes[, c("Genus", "Cohens_d", "Magnitude", "Log2FC")])
    
    p_forest <- ggplot(effect_sizes, aes(x = Cohens_d, y = reorder(Genus, Cohens_d))) +
        geom_vline(xintercept = 0, linetype = "solid", color = "black") +
        geom_point(aes(color = Magnitude), size = 5) +
        geom_segment(aes(xend = 0, yend = Genus, color = Magnitude), linewidth = 1.5) +
        scale_color_manual(values = c("large"="#E63946", "medium"="#F4A261", 
                                       "small"="#457B9D", "negligible"="#ADB5BD")) +
        labs(x = "Cohen's d (HFD vs ND)", y = "",
             title = "Effect Sizes of Five Differential Genera") +
        theme_bw() +
        theme(axis.text.y = element_text(face = "italic", size = 13))
    
    ggsave("effect_sizes_forest_plot.pdf", p_forest, width = 9, height = 4)
    print(p_forest)
    
    write.csv(effect_sizes, "effect_sizes_results.csv", row.names = FALSE)
}

# Save everything
save(five_genera, lag_all_results, effect_sizes, file = "lag_and_effect_size_results.RData")

cat("\nAll analyses complete with fixed names!\n")

# ============================================================================
# RESILIENCE ANALYSIS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("RESILIENCE ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Calculate how quickly communities return to baseline after perturbation
# Using baseline (D1) as reference
baseline_samples <- metadata_clean$Day == 1
baseline_rel <- genus_rel[, baseline_samples]

resilience_data <- data.frame()

for(cage in unique(metadata_clean$cage)) {
    cage_samples <- metadata_clean$cage == cage & metadata_clean$Day > 1
    if(sum(cage_samples) > 0) {
        cage_days <- metadata_clean$Day[cage_samples]
        cage_rel <- genus_rel[, cage_samples, drop = FALSE]
        
        # Get baseline for this cage
        baseline_idx <- which(metadata_clean$cage == cage & metadata_clean$Day == 1)
        if(length(baseline_idx) > 0) {
            baseline <- genus_rel[, baseline_idx]
            
            for(i in 1:ncol(cage_rel)) {
                dist_to_baseline <- vegdist(rbind(baseline, cage_rel[, i]), method = "bray")[1]
                
                resilience_data <- rbind(resilience_data, data.frame(
                    Cage = cage,
                    Diet = unique(metadata_clean$group1[metadata_clean$cage == cage]),
                    Day = cage_days[i],
                    DistanceToBaseline = dist_to_baseline
                ))
            }
        }
    }
}

p_resilience <- ggplot(resilience_data, aes(x = Day, y = DistanceToBaseline, color = Diet)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    labs(x = "Day", y = "Distance to Baseline (D1)",
         title = "Community Resilience: Return to Baseline?",
         subtitle = "Higher distance = less resilient") +
    scale_color_manual(values = diet_colors)

ggsave("resilience_analysis.pdf", p_resilience, width = 8, height = 6)

