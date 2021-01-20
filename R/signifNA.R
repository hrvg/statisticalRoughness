#' Internal function to gracefully handles NA in signif
#' @param value
#' @param digits `numeric` number of significant digits, default to 3
#' @export
#' @keywords internal
signifNA <- function(value, digits = 3){
	if(is.finite(value)){
		return(signif(value, digits))
	} else {
		return("-")
	}
}