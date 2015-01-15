\dontrun{
### Example: ##################################################################
# Model with an un-structured covariance matrix for the random effects of     #
# line and a diagonal residual co-variance matrix.                            #
###############################################################################
rm(list = ls())
library(BGLR)  # only needed to access wheat data. MTM does not depend on BGLR.
data(wheat)
X <- scale(wheat.X, center = TRUE, scale = TRUE)
Y <- wheat.Y
G <- tcrossprod(X)/ncol(X)
nTraits <- ncol(wheat.Y)
df0 <- 3
R2 <- 0.2
S0 <- var(Y)

# Defining the specification of the residual covariance matrix (diagonal).
R0 <- list(type = "DIAG", df0 = rep(df0, nTraits),
           S0 = diag(S0) * (1 - R2) * (df0 + 2))
R0 # Note, when type='DIAG', independent scaled-inv. chi-sq. prior are assigned
# to the diagonal elements, hence, df0 and S0 are vectors, each entry indexes
# one prior.

# Defining the specification of the covariance matrix of random effects
# (un-structured)
KList <- list(list(K = G, COV = list(type = "UN", S0 = diag((1 - R2) * (df0 + 2),
              nTraits), df0 = df0)))
KList[[1]]$COV
# Note for un-structured covariance matrix the prior is Inverse Wishart therefore,
# df0 is a scalar and S0 a pd matrix.

# Fitting the model
fm1 <- MTM(Y = Y, K = KList, resCov = R0, nIter = 1200, burnIn = 200, saveAt = "ex1_")

# Some estimates
fm1$K[[1]]$G  # Estimated genomic covariance matrix
fm1$mu  # Estimated intercepts
fm1$K[[1]]$U  # Predicted genomic values
fm1$resCov$R  # Estimated residual covariance matrix

# Samples and trace plots.  Genomic covariance matrix
G <- read.table("ex1_G_1.dat")
head(G)
plot(G[, 1], type = "o", cex = 0.5, col = 4)  # genomic variance 1st trait.
plot(G[, 2], type = "o", cex = 0.5, col = 4)  # genomic co-variance trait 1 and 2.
# ...
plot(G[, 5], type = "o", cex = 0.5, col = 4)  # genomic variance trait 2.
plot(G[, 6], type = "o", cex = 0.5, col = 4)  # genomic co-variance trait 2 and 3.
library(MCMCpack)
GHat <- xpnd(colMeans(G))

## Residual covariance matrix
R <- read.table("ex1_R.dat")
head(R) # Note, since type='DIAG' several entries (the covariances) are all equal
# to zero.

# Computing Samples for genomic heritabilities (example with trait 1)
H2_1 <- G[, 1]/(G[, 1] + R[, 1])
plot(density(H2_1))
}
