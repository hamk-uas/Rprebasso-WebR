# InitMultiSite nClimID==1 code path corrupts climIDs and breaks ETS computation

## Summary

When all sites share one `climID`, `InitMultiSite` forces `nClimID=2` and reassigns `climIDs[1]=2`. This causes the subsequent ETS loop to fail when called with nSites=1, and produces incorrect climate-to-site mappings with nSites>1.

## Affected code

`R/multiSitePrebas.r`, around line 232 (current master):

```r
if(nClimID == 1){
    nClimID = 2
    climIDs[1] <- 2
    siteInfo[1,2] <- 2
    PAR = matrix(PAR,2,length(PAR),byrow = T)
    TAir = matrix(TAir,2,length(TAir),byrow = T)
    VPD = matrix(VPD,2,length(VPD),byrow = T)
    Precip = matrix(Precip,2,length(Precip),byrow = T)
    CO2 <- matrix(CO2,2,length(CO2),byrow = T)
}
```

Then at line ~286:

```r
for(climID in 1:nClimID){
    nYearsX <- max(nYearsMS[which(siteInfo[,2]==climID)])
    Temp <- TAir[climID,1:(365*nYearsX)]-5
    ...
```

## Observed symptoms

**With nSites=1:** The code sets `climIDs[1]=2`, so `which(siteInfo[,2]==1)` returns no matches. `max()` of an empty vector gives `-Inf`, then `1:(365 * -Inf)` throws:

```
Error in 1:(365 * nYearsX): result would be too long a vector
```

**With nSites=2 (both climID=1):** Site 1 gets its climID changed to 2, so the two sites now point to different weather rows even though they should share the same climate. This produces slightly different results for otherwise identical sites.

## Reproduction

```r
library(Rprebasso)

# Single site — crashes
init <- InitMultiSite(
  nYearsMS = 60L,
  siteInfo = matrix(c(1,1,3,160,0,0,20,3,3,1,0,0), nrow=1),
  multiInitVar = array(c(1,1,7.4,14,0.2,3.1,NA), dim=c(1,7,1)),
  PAR = runif(60*365, 1, 15),
  TAir = runif(60*365, -10, 25),
  VPD = runif(60*365, 0, 2),
  Precip = runif(60*365, 0, 10),
  CO2 = rep(400, 60*365),
  latitude = 60.7
)
# Error in 1:(365 * nYearsX) : result would be too long a vector
```

## Suggested fix

The Fortran `multiPrebas` subroutine and `prebas` should work fine with `nClimID=1`. Remove the artificial inflation and instead handle the single-climate case naturally:

```r
# Remove the if(nClimID == 1) block entirely.
# Ensure weather matrices always have nClimID rows:
if(!is.matrix(PAR)) PAR <- matrix(PAR, nrow=1)
if(!is.matrix(TAir)) TAir <- matrix(TAir, nrow=1)
# ... etc
```

If the Fortran side truly requires `nClimID >= 2`, duplicate the weather row without mutating `climIDs`:

```r
if(nClimID == 1){
    nClimID = 2
    # Do NOT change climIDs — all sites still point to climID row 1
    PAR = rbind(PAR, PAR[1,])
    TAir = rbind(TAir, TAir[1,])
    # ... etc
}
```
