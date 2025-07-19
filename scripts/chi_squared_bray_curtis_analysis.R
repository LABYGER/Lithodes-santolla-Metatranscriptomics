# ========================================================
# Chi-squared Test & Bray-Curtis Dissimilarity in R (Clean Script)
# ========================================================


# --------------------------------------------------------
# 1. Data Preparation
# --------------------------------------------------------
# Example dataset: abundance counts for 2 conditions
abundance_data <- data.frame(
  TaxonA = c(100, 150),
  TaxonB = c(300, 200),
  TaxonC = c(500, 250)
)
rownames(abundance_data) <- c("Condition1", "Condition2")


# --------------------------------------------------------
# 2. Chi-squared Test (Test of Homogeneity)
# --------------------------------------------------------
chisq_test <- chisq.test(abundance_data)

cat("\n=== Chi-squared Test Results ===\n")
print(chisq_test)
cat("\nChi-squared statistic:", chisq_test$statistic, "\n")
cat("p-value:", chisq_test$p.value, "\n")


# --------------------------------------------------------
# 3. Bray-Curtis Dissimilarity (Manual Calculation)
# --------------------------------------------------------
condition1 <- as.numeric(abundance_data["Condition1", ])
condition2 <- as.numeric(abundance_data["Condition2", ])

# Bray-Curtis manual formula
numerator <- sum(abs(condition1 - condition2))
denominator <- sum(condition1 + condition2)
bray_curtis_manual <- numerator / denominator

cat("\n=== Manual Bray-Curtis Dissimilarity ===\n")
cat("Bray-Curtis dissimilarity:", round(bray_curtis_manual, 4), "\n")


# --------------------------------------------------------
# 4. Bray-Curtis Dissimilarity (vegan::vegdist)
# --------------------------------------------------------
library(vegan)
bray_curtis_vegan <- vegdist(abundance_data, method = "bray")

cat("\n=== Bray-Curtis Dissimilarity (vegan::vegdist) ===\n")
print(as.matrix(bray_curtis_vegan))


# --------------------------------------------------------
# 5. Save Results to File
# --------------------------------------------------------
sink("chi_squared_bray_curtis_results.txt")

cat("=== Chi-squared Test Results ===\n")
print(chisq_test)
cat("\nChi-squared statistic:", chisq_test$statistic, "\n")
cat("p-value:", chisq_test$p.value, "\n")

cat("\n=== Manual Bray-Curtis Dissimilarity ===\n")
cat("Bray-Curtis dissimilarity:", round(bray_curtis_manual, 4), "\n")

cat("\n=== Bray-Curtis Dissimilarity (vegan::vegdist) ===\n")
print(as.matrix(bray_curtis_vegan))

sink()
