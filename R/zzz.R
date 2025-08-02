#' @useDynLib PONG
NULL
#' @exportPattern ^[^\\.]
NULL

#' @export kirPredict
NULL
#' @export plot
NULL

#' @export
.onAttach <- function(libname, pkgname) {
  SSE_Flag <- integer(1)

  .C("PONG_Init", SSE_Flag = SSE_Flag, PACKAGE = "PONG")
}

#' @export
.Last.lib <- function(libpath) {
  # Cleanup logic here
}