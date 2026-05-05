library(GWASpoly)

# 1. Load Data
data <- read.GWASpoly(ploidy = 4, 
                      pheno.file = "benchmarks/gwas/new_potato_pheno.csv", 
                      geno.file = "benchmarks/gwas/new_potato_geno.csv", 
                      format = "numeric",
                      n.traits = 1)

# 2. Set Kinship (no LOCO)
data <- set.K(data, LOCO = FALSE)

# 3. Run GWASpoly
data <- GWASpoly(data = data, 
                  models = c("additive", "general"), 
                  traits = "vine.maturity", 
                  params = set.params(MAF = 0.05, geno.freq = 0.95))

# 4. Set threshold using M.eff (default)
data <- set.threshold(data, method = "M.eff", level = 0.05)

# 5. Get QTL
qtl <- get.QTL(data, bp.window = 5e6)
cat("\n============ GWASpoly QTL Results (M.eff, bp.window=5Mb) ============\n")
print(qtl)
cat("\nNumber of QTLs:", nrow(qtl), "\n")

# Also show FDR threshold
data_fdr <- set.threshold(data, method = "FDR", level = 0.05)
qtl_fdr <- get.QTL(data_fdr, bp.window = 5e6)
cat("\n============ GWASpoly QTL Results (FDR, bp.window=5Mb) ============\n")
print(qtl_fdr)
cat("\nNumber of QTLs (FDR):", nrow(qtl_fdr), "\n")

# Show Bonferroni threshold
data_bonf <- set.threshold(data, method = "Bonferroni", level = 0.05)
qtl_bonf <- get.QTL(data_bonf, bp.window = 5e6)
cat("\n============ GWASpoly QTL Results (Bonferroni, bp.window=5Mb) ============\n")
print(qtl_bonf)
cat("\nNumber of QTLs (Bonferroni):", nrow(qtl_bonf), "\n")

# Show ALL hits without windowing  
qtl_all <- get.QTL(data_fdr)
cat("\n============ GWASpoly ALL Hits (FDR, no window) ============\n")
print(qtl_all)
cat("\nNumber of ALL hits (FDR):", nrow(qtl_all), "\n")

# Print the scores for top markers (additive model)
scores <- data@scores$vine.maturity
cat("\n============ Top 20 markers by -log10(p) additive ============\n")
top <- head(scores[order(scores$additive, decreasing=TRUE), ], 20)
print(top)
