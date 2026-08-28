#=============================================================================
# R1-3c & R2-2: VARIANCE PARTITIONING WITH CAGE
#=============================================================================

library(vegan)

# PERMANOVA with cage in model
permanova_full <- adonis2(
    bray_dist ~ group1 * Phase + cage,
    data = metadata_clean,
    permutations = 999,
    by = "terms"
)

print(permanova_full)

# Extract variance components
variance_table <- data.frame(
    Factor = rownames(permanova_full),
    R2 = round(permanova_full$R2, 4),
    P_value = permanova_full$`Pr(>F)`
) %>%
    mutate(Percentage = round(R2 * 100, 1))

print(variance_table)
write.csv(variance_table, "Table_variance_partitioning.csv", row.names = FALSE)

# Plot
library(ggplot2)
p_var <- ggplot(variance_table %>% filter(Factor != "Residual" & Factor != "Total"),
                aes(x = reorder(Factor, R2), y = R2)) +
    geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8, width = 0.6) +
    geom_text(aes(label = paste0(Percentage, "%")), hjust = -0.1, fontface = "bold") +
    coord_flip(ylim = c(0, max(variance_table$R2) * 1.3)) +
    labs(title = "Variance Partitioning of Community Composition",
         x = "", y = "R²") +
    theme_bw()

ggsave("Supplementary_Figure_variance_partitioning.pdf", p_var, width = 8, height = 5)
print(p_var)