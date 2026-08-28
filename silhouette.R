#=============================================================================
# R1-3b: SILHOUETTE ANALYSIS FOR OPTIMAL k
#=============================================================================

library(cluster)
library(ggplot2)

# Calculate silhouette width for k = 2 to 10
sil_df <- data.frame()

for(k in 2:10) {
    set.seed(123)
    pam_fit <- pam(bray_dist, k = k, diss = TRUE)
    sil_df <- rbind(sil_df, data.frame(
        K = k,
        Silhouette = pam_fit$silinfo$avg.width
    ))
}

print(sil_df)

# Plot
p_sil <- ggplot(sil_df, aes(x = K, y = Silhouette)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
    annotate("text", x = 2.3, y = max(sil_df$Silhouette) - 0.02,
             label = "Optimal k = 2", color = "red", hjust = 0, fontface = "bold") +
    labs(title = "Silhouette Analysis for Optimal Cluster Number",
         x = "k", y = "Average Silhouette Width") +
    theme_bw()

ggsave("Supplementary_Figure_silhouette.pdf", p_sil, width = 7, height = 5)
print(p_sil)