## R CMD check results

0 errors | 1 warning | 3 notes

### Warning
- `library(HIBAG)` calls exist inside parallel worker functions in
  kirParallelAttrBagging(). These are required because parallel worker
  nodes start fresh R sessions. HIBAG is listed in Suggests as it is
  an optional enhanced backend.
- `library(parallel)` is called inside kirParallelAttrBagging() only
  when parallel computing is requested by the user.

### Notes
- R_registerRoutines: adding src/init.c conflicts with HIBAG model
  handle lifecycle — dynamic symbol lookup required
- qpdf not available on build system — no PDF vignettes in package
- Startup message: HIBAG attached via tryCatch(attachNamespace())

## Test environments
- Linux x86_64, R 4.4.0

## Downstream dependencies
None — this is a new submission.
