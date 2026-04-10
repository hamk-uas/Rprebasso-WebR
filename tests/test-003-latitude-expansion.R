# Regression test for Issue #003: scalar latitude not expanded to nSites
# Before fix, InitMultiSite() threw:
#   "length of latitude does not match nSites"
# when a single scalar latitude was passed for multiple sites.

library(Rprebasso)

set.seed(42)
nSites <- 2L; nYears <- 5L; ndays <- nYears * 365

initV <- array(NA, dim = c(nSites, 7, 1))
for (i in seq_len(nSites)) initV[i, , 1] <- c(1, 30, 12, 18, 15, 5, NA)

siteInfo <- matrix(c(
  1, 1, 3, 160, 0, 0, 20, 1, 1, 413, 0.45, 0.118, 3,
  2, 2, 3, 160, 0, 0, 20, 1, 1, 413, 0.45, 0.118, 3
), nrow = 2, byrow = TRUE)

weather <- runif(ndays, 1, 15)

# Scalar latitude = 60.7  →  should be replicated to c(60.7, 60.7)
init <- InitMultiSite(
  nYearsMS = rep(nYears, nSites),
  siteInfo = siteInfo,
  multiInitVar = initV,
  PAR    = matrix(rep(weather, each = 2), nrow = 2),
  TAir   = matrix(rep(runif(ndays, -10, 25), each = 2), nrow = 2),
  VPD    = matrix(rep(runif(ndays, 0, 2), each = 2), nrow = 2),
  Precip = matrix(rep(runif(ndays, 0, 10), each = 2), nrow = 2),
  CO2    = matrix(rep(400, ndays * 2), nrow = 2),
  latitude = 60.7  # scalar!
)

stopifnot(
  length(init$latitude) == nSites,
  all(init$latitude == 60.7)
)

res <- multiPrebas(init)
H <- res$multiOut[, , 11, 1, 1]
stopifnot(!any(is.na(H)), all(H > 0))

cat("PASS: test-003-latitude-expansion\n")
