#=============================================================================
# R1-1: TEMPORAL PRECEDENCE - MODEST EFFECT SIZE REPORTING
#=============================================================================

library(dplyr)

# Create effect size summary table
lag_effect_summary <- data.frame(
    Genus = c("Duncaniella", "Taurinivorans", "Lactococcus", "Lepagella", "CAG-873"),
    Peak_Lag = c(-2, 1, -4, 1, 4),
    Peak_r = c(-0.31, 0.31, 0.096, -0.186, 0.053),
    Magnitude = c("small-moderate", "small-moderate", "negligible", "small", "negligible"),
    Interpretation = c(
        "Temporal precedence (abundance precedes weight gain)",
        "Temporal precedence (weight gain precedes abundance)",
        "No clear temporal ordering",
        "Weak temporal association",
        "No clear temporal ordering"
    )
)

print(lag_effect_summary)
write.csv(lag_effect_summary, "Table_lag_effect_sizes.csv", row.names = FALSE)

# Manuscript wording template
cat("\n--- Suggested manuscript wording ---\n")
cat("Cross-correlation analysis revealed temporally ordered associations of modest magnitude:\n")
cat("*Duncaniella* depletion preceded weight gain (peak at lag -2; r = -0.31), while\n")
cat("*Taurinivorans* enrichment followed weight gain (peak at lag +1; r = +0.31).\n")
cat("We interpret this as evidence of temporal precedence, not causation.\n")