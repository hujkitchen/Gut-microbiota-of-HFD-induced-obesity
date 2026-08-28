#=============================================================================
# R2-3: PCoA WITH INDIVIDUAL CAGE TRAJECTORIES
#=============================================================================


library(ggplot2)
library(dplyr)

#=============================================================================
# 1. BUILD FRESH PCOA DATAFRAME FROM AVAILABLE DATA
#=============================================================================

# Get PCoA coordinates if not available
if(!exists("pcoa_df") || !all(c("PC1", "PC2") %in% colnames(pcoa_df))) {
    cat("Rebuilding PCoA dataframe...\n")
    
    # Calculate Bray-Curtis and PCoA
    bray_dist <- vegdist(t(genus_rel), method = "bray")
    pcoa_res <- cmdscale(bray_dist, k = 3, eig = TRUE)
    var_exp <- round(pcoa_res$eig[1:3] / sum(pcoa_res$eig) * 100, 1)
    
    # Build fresh dataframe
    pcoa_df <- data.frame(
        PC1 = pcoa_res$points[,1],
        PC2 = pcoa_res$points[,2],
        PC3 = pcoa_res$points[,3],
        stringsAsFactors = FALSE
    )
}

#=============================================================================
# 2. ADD METADATA TO PCOA DATA
#=============================================================================

# Check sample names match
sample_names <- rownames(pcoa_df)
meta_names <- rownames(metadata_clean)

cat("Sample names match metadata:", all(sample_names %in% meta_names), "\n")

# If no match, try SampleID column
if(!all(sample_names %in% meta_names) && "SampleID" %in% colnames(metadata_clean)) {
    match_idx <- match(sample_names, metadata_clean$SampleID)
} else {
    match_idx <- match(sample_names, meta_names)
}

# Add metadata columns
pcoa_df$Diet <- metadata_clean$group1[match_idx]
pcoa_df$Cage <- metadata_clean$cage[match_idx]
pcoa_df$Day <- as.numeric(metadata_clean$Day[match_idx])
pcoa_df$Phase <- metadata_clean$Phase[match_idx]

# Verify
cat("\nColumns:", paste(colnames(pcoa_df), collapse = ", "), "\n")
cat("Rows:", nrow(pcoa_df), "\n")
cat("Any NA in Diet:", sum(is.na(pcoa_df$Diet)), "\n")
cat("Any NA in Cage:", sum(is.na(pcoa_df$Cage)), "\n")
cat("Any NA in Day:", sum(is.na(pcoa_df$Day)), "\n")

# Remove any rows with NA
pcoa_df <- pcoa_df %>% filter(!is.na(Diet), !is.na(Cage), !is.na(Day))

#=============================================================================
# 3. CREATE TRAJECTORY SEGMENTS
#=============================================================================

# Sort by Cage then Day
pcoa_df <- pcoa_df %>%
    arrange(Cage, Day)

# Calculate next positions for trajectory arrows
pcoa_traj <- pcoa_df %>%
    group_by(Cage, Diet) %>%
    mutate(
        Next_PC1 = lead(PC1),
        Next_PC2 = lead(PC2),
        Next_Day = lead(Day)
    ) %>%
    filter(!is.na(Next_PC1)) %>%
    ungroup()

cat("\nTrajectory segments created:", nrow(pcoa_traj), "\n")

#=============================================================================
# 4. CALCULATE GROUP CENTROIDS
#=============================================================================

centroids <- pcoa_df %>%
    group_by(Diet) %>%
    summarise(
        PC1 = mean(PC1, na.rm = TRUE),
        PC2 = mean(PC2, na.rm = TRUE),
        .groups = 'drop'
    )

#=============================================================================
# 5. PLOT - MAIN FIGURE
#=============================================================================

# X/Y labels
x_lab <- ifelse(exists("var_exp"), paste0("PC1 (", var_exp[1], "%)"), "PC1")
y_lab <- ifelse(exists("var_exp"), paste0("PC2 (", var_exp[2], "%)"), "PC2")

p_trajectory <- ggplot() +
    # Trajectory arrows (light grey)
    geom_segment(
        data = pcoa_traj,
        aes(x = PC1, y = PC2, xend = Next_PC1, yend = Next_PC2, color = Diet),
        alpha = 0.35, linewidth = 0.5,
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        show.legend = FALSE
    ) +
    # Points
    geom_point(
        data = pcoa_df,
        aes(x = PC1, y = PC2, color = Diet, shape = Phase),
        size = 3, alpha = 0.85, stroke = 0.4
    ) +
    # Centroids
    geom_point(
        data = centroids,
        aes(x = PC1, y = PC2, fill = Diet),
        size = 10, shape = 21, color = "black", stroke = 1.8, alpha = 0.9
    ) +
    # Centroid labels
    geom_text(
        data = centroids,
        aes(x = PC1, y = PC2, label = Diet),
        size = 4, fontface = "bold", color = "white"
    ) +
    # Scales
    scale_color_manual(values = c("ND" = "#2E86AB", "HFD" = "#A23B72"), name = "Diet") +
    scale_fill_manual(values = c("ND" = "#2E86AB", "HFD" = "#A23B72"), guide = "none") +
    scale_shape_manual(
        values = c("Acute" = 16, "Early" = 17, "Mid" = 15, "Late" = 18),
        name = "Phase"
    ) +
    # Labels
    labs(
        title = "Gut Microbial Community Trajectories During HFD-Induced Obesity",
        subtitle = paste0(
            "Arrows connect consecutive timepoints within each cage (n = ", 
            length(unique(pcoa_df$Cage)), " cages)\n",
            "Large circles = group centroids | ND (blue) vs HFD (red)"
        ),
        x = x_lab,
        y = y_lab
    ) +
    # Theme
    theme_bw(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10, color = "grey40", lineheight = 1.2),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.minor = element_blank()
    )

ggsave("Revised_Figure1B_trajectories.pdf", p_trajectory, width = 10, height = 8, dpi = 300)
print(p_trajectory)

#=============================================================================
# 6. PLOT - FACETED BY DIET
#=============================================================================

p_traj_facet <- ggplot() +
    # Trajectory arrows
    geom_segment(
        data = pcoa_traj,
        aes(x = PC1, y = PC2, xend = Next_PC1, yend = Next_PC2, color = Diet),
        alpha = 0.4, linewidth = 0.4,
        arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
        show.legend = FALSE
    ) +
    # Points colored by Day
    geom_point(
        data = pcoa_df,
        aes(x = PC1, y = PC2, fill = Day),
        shape = 21, size = 3, color = "white", stroke = 0.3, alpha = 0.9
    ) +
    # Facet by Diet
    facet_wrap(~ Diet, ncol = 2) +
    # Scales
    scale_fill_viridis_c(option = "plasma", name = "Day") +
    # Labels
    labs(
        title = "Community Trajectories by Diet Group",
        subtitle = "Arrows connect consecutive timepoints within each cage",
        x = x_lab,
        y = y_lab
    ) +
    # Theme
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom"
    )

ggsave("Revised_Figure1B_trajectories_faceted.pdf", p_traj_facet, 
       width = 12, height = 6, dpi = 300)
print(p_traj_facet)

#=============================================================================
# 7. SAVE
#=============================================================================

save(p_trajectory, p_traj_facet, pcoa_df, pcoa_traj, centroids,
     file = "revised_figure1B_data.RData")

cat("\nDone! Files saved:\n")
cat("  - Revised_Figure1B_trajectories.pdf\n")
cat("  - Revised_Figure1B_trajectories_faceted.pdf\n")