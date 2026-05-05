# 1. Load GWASpoly
library(GWASpoly)

# 2. Load Data
geno_file <- "benchmarks/gwas/new_potato_geno.csv"
pheno_file <- "benchmarks/gwas/new_potato_pheno.csv"

# GWASpoly 2.x read.GWASpoly
# max.missing is 0.1 by default in read.GWASpoly if not specified, 
# but let's check the function signature if possible.
# For now, let's just use defaults.
data <- read.GWASpoly(ploidy = 4, 
                      pheno.file = pheno_file, 
                      geno.file = geno_file, 
                      format = "numeric",
                      n.traits = 1)

# 3. Set Kinship
# 4. Set Kinship
data <- set.K(data, LOCO = FALSE)

# Estimate variance components for the trait
# We can use GWASpoly's internal function or just run a dummy scan
print("Estimating variance components manually...")
params <- set.params(MAF=0.05, geno.freq=0.95, fixed="env", fixed.type="factor")
# For GWASpoly 2.x, the variance components are usually estimated during GWAS()
# Let's run a quick scan on 10 markers to get them
data_quick <- GWASpoly(data = data, 
                  models = "additive", 
                  traits = "vine.maturity", 
                  params = params,
                  n.core = 1)

print("Slots in data_quick:")
print(slotNames(data_quick))

# In GWASpoly 2.x, the results are in @scores
# Variance components might be in @params[[trait]]
print("Params for vine.maturity:")
print(data_quick@params$vine.maturity)

# If not there, check if it's in the attributes of the scores
print("Attributes of scores:")
print(attributes(data_quick@scores$vine.maturity))


# Inspect Kinship Matrix
print("Kinship Matrix (first 5x5):")
if (is.list(data_quick@K)) {
  K_mat <- data_quick@K[[1]]
} else {
  K_mat <- data_quick@K
}
print(K_mat[1:5, 1:5])
print(paste("Mean Diagonal of K:", mean(diag(K_mat))))

# Quit early
quit(save="no")


# Check fixed effects
print("Fixed effects used in the model:")
print(data@fixed)


# Print slots to debug
print("Slots available in the data object:")
print(slotNames(data))

# Print Variance Components (GWASpoly 2.x style)
# Usually they are in @Vg and @Ve but let's be safe
if ("Vg" %in% slotNames(data)) {
  print("Variance Components:")
  print(data@Vg)
  print(data@Ve)
} else {
  print("Vg slot not found. Printing the whole object structure for the first trait:")
  # Sometimes it's in a list
}


# 5. Get Results
results <- data@scores$vine.maturity

# 6. Threshold
n_markers <- nrow(data@geno)
threshold <- -log10(0.05 / n_markers)
print(paste("Bonferroni threshold (-log10p):", threshold))

# 7. Print Top Results
sig_hits_mask <- rowSums(results > threshold, na.rm=TRUE) > 0
if (any(sig_hits_mask)) {
  print("Significant Hits (> Bonferroni):")
  print(head(results[sig_hits_mask, ], 20))
} else {
  print("No significant hits found above Bonferroni.")
}

print("Top 10 Markers by -log10(p) (Additive Model):")
print(head(results[order(results$additive, decreasing = TRUE), ], 10))
