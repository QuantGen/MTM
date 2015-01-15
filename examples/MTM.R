\dontrun{
library(BGLR)

data(wheat)
G <- tcrossprod(scale(wheat.X, center = TRUE, scale = TRUE)) / ncol(wheat.X)

nTraits <- ncol(wheat.Y)
df0 <- 3
R2 <- 0.5

# K is similar to ETA in BGLR, but only for RKHS in the case of this program.
# One can provide K or EVD via V and d, I believe. Use much longer chains, this
# is just an example!
fm <- MTM(Y = wheat.Y,
          K = list(list(K = G, COV = list(type = "UN", df0 = (df0 + nTraits),
                S0 = diag(R2 * (df0 + 2 * nTraits + 1), nTraits)))),
          resCov = list(type = "DIAG", df0 = rep(df0, nTraits),
                S0 = diag(R2 * df0/(df0 + 2), nTraits)),
          nIter = 500,
          burnIn = 100)

# Some estimates
fm$K[[1]]$G  # Estimated genomic covariance matrix
fm$K[[1]]$U  # Predicted genomic values
fm$resCov$R  # Estimated residual covariance matrix
}
