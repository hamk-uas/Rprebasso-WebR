# Regression test for Issue #002: nClimID==1 crash
# Before fix, InitMultiSite() with a single climate ID and nSites==1
# crashed because it mutated climIDs[1] to 2 and then tried to index
# weather data that didn't exist for climID 2.

library(Rprebasso)

set.seed(42)
nYears <- 5L; ndays <- nYears * 365

# --- Case A: single site, single climID (the crashing case) ---
initV1 <- array(NA, dim = c(1, 7, 1))
initV1[1, , 1] <- c(1, 30, 12, 18, 15, 5, NA)

siteInfo1 <- matrix(c(1, 1, 3, 160, 0, 0, 20, 1, 1, 413, 0.45, 0.118, 3),
                    nrow = 1)

init1 <- InitMultiSite(
  nYearsMS   = nYears,
  siteInfo   = siteInfo1,
  multiInitVar = initV1,
  PAR    = matrix(runif(ndays, 1, 15), nrow = 1),
  TAir   = matrix(runif(ndays, -10, 25), nrow = 1),
  VPD    = matrix(runif(ndays, 0, 2), nrow = 1),
  Precip = matrix(runif(ndays, 0, 10), nrow = 1),
  CO2    = matrix(rep(400, ndays), nrow = 1),
  latitude = 60.7
)

# climIDs must NOT be mutated
stopifnot(all(init1$siteInfo[, 2] == 1))

res1 <- multiPrebas(init1)
H1 <- res1$multiOut[1, , 11, 1, 1]
stopifnot(!any(is.na(H1)), all(H1 > 0))

cat("PASS: test-002a (single site, single climID)\n")

# --- Case B: two sites sharing one climID (was silently corrupted) ---
initV2 <- array(NA, dim = c(2, 7, 1))
initV2[1, , 1] <- c(1, 30, 12, 18, 15, 5, NA)
initV2[2, , 1] <- c(1, 30, 12, 18, 15, 5, NA)

siteInfo2 <- matrix(c(
  1, 1, 3, 160, 0, 0, 20, 1, 1, 413, 0.45, 0.118, 3,
  2, 1, 3, 160, 0, 0, 20, 1, 1, 413, 0.45, 0.118, 3
), nrow = 2, byrow = TRUE)

init2 <- InitMultiSite(
  nYearsMS   = c(nYears, nYears),
  siteInfo   = siteInfo2,
  multiInitVar = initV2,
  PAR    = matrix(runif(ndays, 1, 15), nrow = 1),
  TAir   = matrix(runif(ndays, -10, 25), nrow = 1),
  VPD    = matrix(runif(ndays, 0, 2), nrow = 1),
  Precip = matrix(runif(ndays, 0, 10), nrow = 1),
  CO2    = matrix(rep(400, ndays), nrow = 1),
  latitude = c(60.7, 60.7)
)

stopifnot(all(init2$siteInfo[, 2] == 1))

res2 <- multiPrebas(init2)
HA <- res2$multiOut[1, , 11, 1, 1]
HB <- res2$multiOut[2, , 11, 1, 1]
stopifnot(!any(is.na(HA)), all(HA > 0))
# Same init + same weather → same result
stopifnot(all(abs(HA - HB) < 0.01))

cat("PASS: test-002b (two sites, shared climID, identical results)\n")
