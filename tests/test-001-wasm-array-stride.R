# Regression test for Issue #001: non-contiguous array strides in Fortran
# Before fix, multiPrebas passed non-contiguous array sections directly to
# prebas(), which receives explicit-shape dummy arrays. Some compilers create
# contiguous temporaries for this call pattern, but the WebAssembly/flang build
# did not, so key state variables were read with the wrong stride.
#
# Native R cannot reproduce the compiler-specific flang failure directly, so
# this regression checks the variables reported in the issue (H, D, GPPtrees,
# NEP) and verifies that identical sites stay numerically aligned.

library(Rprebasso)

set.seed(42)
nSites <- 3L; nYears <- 5L; ndays <- nYears * 365

initV <- array(NA, dim = c(nSites, 7, 1))
for (i in seq_len(nSites)) initV[i, , 1] <- c(1, 30, 12, 18, 15, 5, NA)

siteInfo <- matrix(nrow = nSites, ncol = 13)
for (i in seq_len(nSites)) {
  siteInfo[i, ] <- c(i, i, 3, 160, 0, 0, 20, 1, 1, 413, 0.45, 0.118, 3)
}

PAR_d  <- runif(ndays, 1, 15)
TAir_d <- runif(ndays, -10, 25)
VPD_d  <- runif(ndays, 0, 2)
Prec_d <- runif(ndays, 0, 10)

init <- InitMultiSite(
  nYearsMS = rep(nYears, nSites),
  siteInfo = siteInfo,
  multiInitVar = initV,
  PAR    = matrix(rep(PAR_d, each = nSites), nrow = nSites),
  TAir   = matrix(rep(TAir_d, each = nSites), nrow = nSites),
  VPD    = matrix(rep(VPD_d, each = nSites), nrow = nSites),
  Precip = matrix(rep(Prec_d, each = nSites), nrow = nSites),
  CO2    = matrix(rep(400, ndays * nSites), nrow = nSites),
  latitude = rep(60.7, nSites)
)

res <- multiPrebas(init)

H <- res$multiOut[, , 11, 1, 1]
D <- res$multiOut[, , 12, 1, 1]
BA <- res$multiOut[, , 13, 1, 1]
GPPtrees <- res$multiOut[, , 44, 1, 1]
NEP <- res$multiOut[, , 46, 1, 1]

stopifnot(!any(is.na(H)), all(H > 0))
stopifnot(!any(is.na(D)), all(D > 0))
stopifnot(!any(is.na(BA)), all(BA > 0))
stopifnot(!any(is.na(GPPtrees)), all(GPPtrees > 0))
stopifnot(!any(is.na(NEP)))

# Identical sites and identical weather should produce matching site trajectories.
stopifnot(all(abs(H[1, ] - H[2, ]) < 1e-8))
stopifnot(all(abs(H[1, ] - H[3, ]) < 1e-8))
stopifnot(all(abs(D[1, ] - D[2, ]) < 1e-8))
stopifnot(all(abs(BA[1, ] - BA[2, ]) < 1e-8))
stopifnot(all(abs(GPPtrees[1, ] - GPPtrees[2, ]) < 1e-8))
stopifnot(all(abs(NEP[1, ] - NEP[2, ]) < 1e-8))

# The issue report called out H~=1 and NEP=NaN under Wasm corruption.
stopifnot(all(H[1, ] > 10))
stopifnot(all(is.finite(NEP[1, ])))

cat("PASS: test-001-wasm-array-stride\n")
