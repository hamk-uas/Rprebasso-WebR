# C_multiPrebas passes non-contiguous Fortran array sections to prebas — breaks under WASM

## Summary

`C_multiPrebas.f90` passes non-contiguous array sections (e.g. `initVar(i,:,1:nLayers(i))`) to the `prebas()` subroutine. Under native gfortran this works because the compiler generates temporary contiguous copies. Under WASM (flang / Emscripten), no temporary is created for these calls, so `prebas()` reads flat memory with the wrong stride and gets scrambled values.

The important nuance is that this is an explicit-shape dummy argument problem, not the assumed-shape case described by Flang's array repacking flags. `prebas()` declares `initVar(7, nLayers)`, `siteInfo(11)`, `output(nYears, nVar, nLayers, 2)`, etc., and `initBiomasses()` declares `initVar(7)`. Because of that, `-frepack-arrays` is not a reliable fix for this issue as the code currently stands.

## Affected code

`src/C_multiPrebas.f90`, line ~113:

```fortran
call prebas(..., initVar(i,:,1:nLayers(i)), ...)
```

`initVar` has shape `(nSites, 7, maxNlayers)`. When `i` fixes the first dimension, the resulting section `initVar(i,:,1:nLayers(i))` has stride `nSites` in memory (column-major). The `prebas` subroutine declares `initVar(7, nLayers)` and reads it as if contiguous.

Other per-site extractions in the same call have the same problem: `siteInfo(i,:)`, `output(...)`, etc.

Relevant callee declarations:

- `src/B_prebas.f90`: `initVar(7, nLayers)`, `siteInfo(11)`, `output(nYears, nVar, nLayers, 2)`
- `src/A_routines.f90`: `initVar(7)` in `initBiomasses()`

## Observed symptoms

With nSites=2, `prebas()` reads initVar as contiguous flat memory starting at the wrong offset:

| Variable | Expected | Actual (WASM) | Explanation |
|----------|----------|---------------|-------------|
| H (initVar(3,1)) | 7.4495 | 1 | Reads flat[3] = age value |
| D (initVar(4,1)) | 14.348 | 1 | Reads flat[4] = age of site 2 |
| GPPtrees (var 44) | ~3.0 | 0 | Growth fails with H=1 |
| NEP (var 46) | ~0.4 | NaN | Derived from broken GPP/Rh |

With nSites=1 there is no stride issue (stride=1), so the same code works.

## Verified with native R

```r
# nSites=2, identical sites
result <- multiPrebas(init)
out <- result$multiOut
out[1, 1:5, 11, 1, 1]  # H site1 yr1-5: 7.52 7.59 7.66 7.72 7.78  ← correct
out[1, 1:5, 46, 1, 1]  # NEP site1 yr1-5: 0.42 0.41 0.39 0.38 0.38  ← correct
```

Same input in WebR gives H=1, NEP=NaN.

## Compiler flags

Flang does document shipped repacking-related flags in its current command-line reference, including:

- `-frepack-arrays`
- `-fstack-repack-arrays`
- `-frepack-arrays-contiguity=whole|innermost`

However, those flags are documented for non-contiguous assumed-shape dummy arrays. They do not provide a robust compiler-option-only fix for this issue, because the failing callees here use explicit-shape dummy arrays.

Two practical consequences follow:

1. `-frepack-arrays` is not something we should rely on to fix issue 001 in the current API shape.
2. `-fstack-repack-arrays` would be the wrong direction for WebAssembly anyway, because stack-backed temporaries are more fragile in the WebR/Wasm environment than heap-backed ones.

## Suggested fix

In `C_multiPrebas.f90`, explicitly copy non-contiguous sections into local contiguous temporaries before passing to `prebas()`:

```fortran
! Before the prebas call:
real(kind=8) :: initVar_tmp(7, maxNlayers)
initVar_tmp(:, 1:nLayers(i)) = initVar(i, :, 1:nLayers(i))

call prebas(..., initVar_tmp, ...)
```

Apply the same pattern to all array arguments that use per-site slicing.

Implementation note: large temporaries should be heap-backed (`allocatable`) rather than large automatic stack arrays. That keeps the fix safe for WebAssembly, where stack size is typically much tighter than on native builds.

This is the approach now implemented in `src/C_multiPrebas.f90`: explicit copy-in/copy-out with heap-backed work buffers for the large per-site arrays.

Alternatively, restructure the data layout so that `nSites` is the last dimension rather than the first, making per-site slices naturally contiguous.
