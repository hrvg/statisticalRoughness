#' This wrapper function derives the one-dimensional `HHCF()` for all columns or rows of a matrix
#' @param mat a `matrix`
#' @param dr `numeric`, spacing of the values along the axis
#' @param margin indicate the direction, 1 = row, 2 = column, passed to `apply`'s `MARGIN`
#' @param detrend `logical` indicates if the data need to be detrended, default to `TRUE`
#' @param get_autocorr `logical`, indicaites if the auto-correlation length is calculated
#' @export
#' @keywords zeta
get_hhcf <- function(mat, dr, margin = 1, detrend = TRUE, get_autocorr = FALSE, limlen = 30){
	hhcf <- plyr::alply(mat, .margins = margin, .fun = function(row){
		ind <- which(!is.na(row))
		if (length(ind) > limlen) {
			if (detrend){
				x <- seq_along(row) * dr
				len <- length(row)
				fit <- lm(row ~ x)
				row <- fit$residuals
			} else {
				row <- row[ind]
			}
			return(HHCF(row))
		} else {
			return(NA)
		}
	})

	if (get_autocorr){
		autocorr_len <- plyr::alply(mat, .margins = margin, .fun = function(row){
			ind <- which(!is.na(row))
			if (length(ind) > limlen) {
				if (detrend){
					x <- seq_along(row) * dr
					len <- length(row)
					fit <- lm(row ~ x)
					row <- fit$residuals
				} else {
					row <- row[ind]
				}
				ACF <- stats::acf(row, plot = FALSE, type = "correlation", lag.max = length(row) - 1)
				return(ACF$lag[utils::head(which(ACF$acf <= 1 / exp(1)), 1)])
			} else {
				return(NA)
			}
		})
		autocorr_len <- unname(unlist(autocorr_len)) * dr
	} else {
		autocorr_len <- NULL
	}

	max_length <- max(sapply(hhcf, length))
	hhcf <- lapply(hhcf, function(row){
		length(row) <- max_length
		return(row)
	})
	hhcf <- do.call(rbind, hhcf)
	hhcf <- data.frame(hhcf)
	# hhcf <- hhcf[rowSums(is.na(hhcf)) != ncol(hhcf), ]
	len <- apply(hhcf, MARGIN = 1, FUN = function(row) length(stats::na.omit(row)))
	w <- apply(mat, MARGIN = margin, FUN = function(row){
		ind <- which(!is.na(row))
		if (length(ind) > limlen) {
			return(RMS_roughness(row, dr))
		} else{
			return(NA)
		}
	})
	return(list(hhcf = hhcf, len = len, autocorr_len = autocorr_len, rms = w))
}

#' Internal; derives a one-dimensional height-height correlation function
#' @param row array of values
#' @return array with the values of the height-height correlation function
#' @export
#' @keywords zeta
HHCF <- function(row){
	len <- length(row) - 1
	dr <- seq(1, len)
	diffs <- lapply(dr, function(r) diff(row, lag = r))
	hhcf <- sapply(diffs, function(d) sqrt(sum(d^2) / length(d)))
	return(hhcf)
}

RMS_roughness <- function(row, dr){
	dr <- seq_along(row) * dr
	len <- length(row)
	fit <- lm(row ~ dr)
	row <- fit$residuals
	avg <- sum(row, na.rm = FALSE)  / len
	w <- sum((row - avg)^2, na.rm = FALSE) / (len - 1) # this is faster than mean()
	return(sqrt(w))
}


#' This wrapper function derives the one-dimensional `HHCF()` for all columns or rows of a matrix
#' @param mat a `matrix`
#' @param dr `numeric`, spacing of the values along the axis
#' @param margin indicate the direction, 1 = row, 2 = column, passed to `apply`'s `MARGIN`
#' @param limlen `numeric` the minimum length of data
#' @export
#' @keywords zeta
get_hhcf_ <- function(mat, dr, margin = 1, limlen = 30){
	hhcf <- plyr::alply(mat, .margins = margin, .fun = function(row){
		ind <- which(!is.na(row))
		if (length(ind) > limlen) {
			x <- seq_along(row) * dr
			len <- length(row)
			fit <- lm(row ~ x)
			row <- row - fit$fitted.values
			ACV <- stats::acf(row, plot = FALSE, type = "covariance", demean = FALSE, lag.max = length(row) - 1, na.action = na.pass)
			W <- ACV$acf[1]
			HHCF <- sqrt(2 * W - 2 * ACV$acf)[-1]
			ACF <- ACV$acf / W
			xi <- ACV$lag[utils::head(which(ACF <= 1 / exp(1)), 1)]
			return(list(hhcf = HHCF, w = sqrt(W), xi = xi))
		} else {
			return(NA)
		}
	})
	ind <- which(is.na(hhcf))
	if (length(ind) > 1) {
		hhcf <- rlist::list.remove(hhcf, ind)
	}
	w <- sapply(hhcf, function(elmt) elmt[[2]])
	xi <- sapply(hhcf, function(elmt) elmt[[3]]) * dr
	hhcf <- lapply(hhcf, function(elmt) elmt[[1]])
	hhcf <- do.call(rbind, hhcf) # one hhcf per line
	hhcf <- data.frame(hhcf)
	return(list(hhcf = hhcf, autocorr_len = xi, rms = w))
}
