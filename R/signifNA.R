#' Internal function to gracefully handles NA in signif
#' @param value value to be transformed
#' @param digits `numeric` number of significant digits, default to 3
#' @export
#' @keywords internal
signifNA <- function(value, digits = 3) {
  if (is.finite(value)) {
    return(signif(value, digits))
  } else {
    return("-")
  }
}
