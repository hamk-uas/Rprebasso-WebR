# InitMultiSite does not expand scalar `latitude` to nSites length

## Summary

When a scalar `latitude` value is passed to `InitMultiSite`, it is stored as-is (length 1) rather than being expanded to length `nSites`. Downstream code that indexes `latitude[s]` for `s > 1` gets `NA`, which causes the `.Fortran("multiPrebas",...)` call to fail with:

```
Error in multiPrebas(d): NA/NaN/Inf in foreign function call (arg 54)
```

Argument 54 in the `.Fortran` call is `latitude = as.double(multiSiteInit$latitude)`.

## Affected code

`R/multiSitePrebas.r`, around line 180:

```r
if(all(is.na(latitude))){
    warning("latitude was not provided. a default value of 62 was used.")
    latitude = rep(62, nSites)
}
```

When `latitude` is not `NA`, no expansion happens. Compare with other scalar parameters that are correctly expanded:

```r
if (length(defaultThin) == 1) defaultThin = as.double(rep(defaultThin, nSites))
if (length(ClCut) == 1) ClCut = as.double(rep(ClCut, nSites))
if (length(energyCut) == 1) energyCut = as.double(rep(energyCut, nSites))
```

## Suggested fix

```r
if(all(is.na(latitude))){
    warning("latitude was not provided. a default value of 62 was used.")
    latitude = rep(62, nSites)
} else if (length(latitude) == 1) {
    latitude = rep(latitude, nSites)
}
```
