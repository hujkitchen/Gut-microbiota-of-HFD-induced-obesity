#=============================================================================
# R2-5: BIAS-CORRECTED LOG2FC TRAJECTORIES
#=============================================================================

library(ggplot2)
library(dplyr)

# Calculate per-day Log2FC for the five genera
day_fc_data <- data.frame()

for(fg in five_genera) {
    # Get abundance
    if(fg %in% rownames(genus_rel)) {
        abund <- genus_rel[fg, ]
    } else {
        abund <- genus_rel[paste0("g__", fg), ]
    }
    
    for(d in sort(unique(metadata_clean$Day))) {
        day_idx <- which(metadata_clean$Day == d)
        nd_vals <- abund[day_idx[metadata_clean$group1[day_idx] == "ND"]]
        hfd_vals <- abund[day_idx[metadata_clean$group1[day_idx] == "HFD"]]
        
        if(length(nd_vals) > 0 && length(hfd_vals) > 0) {
            log2fc <- log2((mean(hfd_vals) + 1e-10) / (mean(nd_vals) + 1e-10))
            
            day_fc_data <- rbind(day_fc_data, data.frame(
                Genus = fg, Day = d, Log2FC = log2fc
            ))
        }
    }
}

# Plot
p_fc <- ggplot(day_fc_data, aes(x = Day, y = Log2FC, color = Genus)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Lactococcus" = "#E41A1C",
                                   "Duncaniella" = "#377EB8",
                                   "Taurinivorans" = "#4DAF4A",
                                   "Lepagella" = "#984EA3",
                                   "CAG-873" = "#FF7F00")) +
    scale_x_continuous(breaks = sort(unique(day_fc_data$Day))) +
    labs(title = "Temporal Log2FC (HFD vs ND) for Five Focal Genera",
         subtitle = "Positive = enriched in HFD; negative = depleted in HFD",
         x = "Day", y = "Log2 Fold Change") +
    theme_bw() +
    theme(legend.text = element_text(face = "italic"))

ggsave("Revised_Figure2C_log2fc.pdf", p_fc, width = 10, height = 6)
print(p_fc)

# Save
write.csv(day_fc_data, "Table_bias_corrected_log2fc.csv", row.names = FALSE)