library(GWASpoly)
library(rrBLUP)

# Load data
pheno <- read.csv("benchmarks/gwas/new_potato_pheno.csv", check.names=FALSE)
cat("Pheno Cols:", colnames(pheno), "\n")
geno <- read.csv("benchmarks/gwas/new_potato_geno.csv", row.names=1, check.names=FALSE)

# Basic setup to get K
# We need to filter exactly like BioCenicana/GWASpoly
# MAF 0.05, max.missing 0.1, max.geno.freq 0.95
# But let's just get the Kinship from ZZ'

# GWASpoly scaling: K = ZZ' / mean(diag(ZZ'))
# dosages (0..4)
M <- as.matrix(geno[,-c(1,2)]) # Assuming first 2 are CHR, POS? No, geno in new_potato_geno.csv has Marker as rowname
# Let's check structure
# head(geno)
#        Marker Chrom Position Sample1 Sample2 ...
# Marker1 ...

M <- as.matrix(geno[,-c(1,2)])
Z <- t(M) # Samples as rows
# Impute missing with mean
Z[is.na(Z)] <- mean(Z, na.rm=TRUE)
Z <- scale(Z, center=TRUE, scale=FALSE)

ZZt <- Z %*% t(Z)
K <- ZZt / mean(diag(ZZt))

# Match samples and handle NAs
common <- intersect(pheno$id, rownames(Z))
cat("Common samples:", length(common), "\n")
pheno_sub <- pheno[match(common, pheno$id), ]
keep <- !is.na(pheno_sub$vine.maturity)
y <- pheno_sub$vine.maturity[keep]
cat("Y size:", length(y), "\n")
K_sub <- K[common[keep], common[keep]]
X <- matrix(1, length(y), 1)

target_samples <- c("AF5392-8", "AF5393-1", "AF5445-2")
cat("K_sub specific samples block:\n")
print(K_sub[target_samples, target_samples])

# Solve mixed model
out <- mixed.solve(y=y, X=X, K=K_sub)
cat("GWASpoly REML Results:\n")
cat("VarG:", out$Vu, "\n")
cat("VarE:", out$Ve, "\n")
cat("Delta (Ve/Vg):", out$Ve / out$Vu, "\n")
