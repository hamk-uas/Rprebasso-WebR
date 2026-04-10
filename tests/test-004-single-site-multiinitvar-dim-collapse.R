library(Rprebasso)

set.seed(42)

nYears <- 60L
nDays <- nYears * 365

siteInfo <- matrix(c(1, 1, 3, 160, 0, 0, 20, 3, 3, 413, 0.45, 0.118, 3),
                   nrow = 1)

multiInitVar <- array(0, dim = c(1, 7, 3))
multiInitVar[1, 1, ] <- c(1, 2, 3)
multiInitVar[1, 2, ] <- c(1, 1, 1)
multiInitVar[1, 3, ] <- c(7.4495, 5.7, 4.3)
multiInitVar[1, 4, ] <- c(14.348, 10, 9)
multiInitVar[1, 5, ] <- c(0.2, 0.1, 0.0)
multiInitVar[1, 6, ] <- c(3.1, 2.3, 2.3)
multiInitVar[1, 7, ] <- c(NA, NA, NA)

PAR <- matrix(runif(nDays, 1, 15), nrow = 1)
TAir <- matrix(runif(nDays, -10, 25), nrow = 1)
VPD <- matrix(runif(nDays, 0, 2), nrow = 1)
Precip <- matrix(runif(nDays, 0, 10), nrow = 1)
CO2 <- matrix(rep(400, nDays), nrow = 1)

init <- InitMultiSite(
  nYearsMS = nYears,
  siteInfo = siteInfo,
  multiInitVar = multiInitVar,
  PAR = PAR,
  TAir = TAir,
  VPD = VPD,
  Precip = Precip,
  CO2 = CO2,
  latitude = 60.7,
  pCROBAS = pCROB,
  pPRELES = pPREL,
  defaultThin = 0,
  ClCut = 0,
  inDclct = NA,
  inAclct = NA
)

stopifnot(
  length(dim(init$multiInitVar)) == 3,
  identical(dim(init$multiInitVar), c(1L, 7L, 3L)),
  all(init$multiInitVar[1, 2, ] >= 1)
)

res <- multiPrebas(init)
H <- res$multiOut[1, , 11, , 1]

stopifnot(
  !any(is.na(H[, 1:2, drop = FALSE])),
  all(H[, 1:2, drop = FALSE] > 0),
  all(is.finite(H[, 3]))
)

cat("PASS: test-004-single-site-multiinitvar-dim-collapse\n")