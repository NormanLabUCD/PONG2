## Resubmission v1.0.1

This is a resubmission addressing all CRAN reviewer comments:

1. Added references in DESCRIPTION: Zheng et al. (2014) <doi:10.1016/j.ajhg.2013.12.015>
2. Added \value tags to all .Rd files including PONG2-internal.Rd
3. Added executable examples to kirPredict.Rd and kirParallelAttrBagging.Rd using built-in PONG2_example data
4. Added PONG2_example dataset documentation in man/PONG2_example.Rd
5. Added PONG2-R-api.Rmd vignette with fully executable R code (eval=TRUE)
6. Replaced cat()/print() with message() throughout R/PONG.r
7. Replaced writing to package install directory with tools::R_user_dir() in .get_model_path()
8. Removed setup_PONG2() which wrote to ~/.local/bin/ without consent
9. Added Xiuwen Zheng as contributor (ctb, cph) in Authors@R
10. Added tools to Imports for tools::R_user_dir()
11. Added importFrom("utils", "capture.output") to NAMESPACE

## R CMD check results
- 0 errors
- 1 warning: 'qpdf' is needed for checks on size reduction of PDFs (environment issue only)
- 1 note: unable to verify current time (network issue on cluster)
