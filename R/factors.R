#' Get factor of an integer
#' @param x, numeric, converted to integer
#' @return list of factors
#' @export
#' @keywords factors
get_factors <- function(x) {
    x <- as.integer(x)
    div <- seq_len(abs(x))
    factors <- div[x %% div == 0L]
    return(factors)
}

#' Get a list of upper spatial scale at which the roughness coefficient is calculated
#' @param Lmax maximum value of L
#' @param factor_limit a threshold of the minimum length of factors for a scale to be accepted
#' @param only default to `NULL`, if `only == 610` only factors of 6 and 10 are considered for a quick search, if `only == 12` only factors of 12 are considered for the fastest search
#' @param logfilter logical, default to `TRUE`, if `TRUE` the returned scale are logaritmically spaced
#' @param len integer, number of scales returned, default to 32, len only works if `logfilter == TRUE`
#' @return a list with two elements `allRs` and `allL`
#' @export
#' @keywords factors
get_all_R_L <- function(Lmax, factor_limit, only = NULL, logfilter = TRUE, len = 32){
	allRs  <- list()
	allL  <- c()
	Ls <- seq(Lmax)
	if (!is.null(only)){
		if (only == 610){
			i610 <- which(Ls %% 6 == 0 | Ls %% 10 == 0)
			Ls <- Ls[i610]
		}
		if (only == 12){
			i12 <- which(Ls %% 12 == 0)
			Ls <- Ls[i12]
		}
	}
	for (L in Ls){
		Rs <- get_factors(L)
		Rs <- utils::tail(Rs, -1)
		Rs <- utils::head(Rs, -1)
		if (length(Rs) >= factor_limit){
			allRs <- append(allRs, list(Rs))
			allL <- append(allL, L)
		}
	}
	if (logfilter){
		if (length(allL) < len) stop('Number of factors return is lower than requested.')
		ind <- 1
		k <- 0
		while (length(ind) != len){
			ind <- unique(as.integer(pracma::logseq(1, length(allL), n = len + k)))
			k <- k + 1
		}	
		return(list(allRs = allRs[ind], allL = allL[ind]))
	} else {
		return(list(allRs = allRs, allL = allL))
	}

}