#' Wrapper function around the computation of Hurst coefficients from box-counting and log-log regression
#' @param dem initial DEM
#' @param L upper scale region size
#' @param R list of scales (aggregation factors)
#' @return a RasterStack with four layers q, H, p-value and R squared values
#' @export
compute_Hurst_rasters <- function(dem, L, R){
	L_raster <- terra::aggregate(dem, fact = L, fun = mean, na.rm = TRUE)
	mean_res <- get_mean_res(L_raster, L)
	mean_sd_rasters <- get_mean_sd_rasters(dem, L, R)
	Hurst_rasters <- compute_Hurst_rasters_internal(mean_sd_rasters, R, mean_res, L_raster)
	return(Hurst_rasters)
}

#' Compute Hurst rasters from a log-log regression between scales and mean standard deviation
#' @param mean_sd_rasters list of rasters of mean standard deviation computed at R scales
#' @param R list of scales (aggregation factors)
#' @param mean_res numeric, resolution to compute scales in meters
#' @param L_raster target raster at scale L
#' @param min_res minimum resolution in elevation, 1.87 m from Haneberg, 2006
#' @return a RasterStack with four layers q, H, p-value and R squared values
#' @import foreach
#' @import doFuture
#' @importFrom stats lm lm.fit pt
#' @export
compute_Hurst_rasters_internal <- function(mean_sd_rasters, R, mean_res, L_raster, min_res = 1.87){
	x <- logsd <- r <- id <- variable <- NULL
	raster_values <- foreach(x = mean_sd_rasters, r = R, .combine = rbind) %do% {
		data.frame(
			id = seq(terra::ncell(x)),
			logsd = unname(terra::values(x))
		) %>% 
		dplyr::mutate(logsd = dplyr::case_when(logsd < min_res ~ min_res, TRUE ~ logsd)) %>%
		dplyr::mutate(logsd = log10(logsd)) %>%
		dplyr::mutate(logscale = log10(r * mean_res))
	}

	.dofastlm <- function(df){
		x <- cbind(1, df$logscale)
		y <- df$logsd
		if(all(is.finite(y))){
			z <- lm.fit(x, y)
			p <- z$rank
			p1 <- 1L:p
	   		R <- chol2inv(z$qr$qr[p1, p1, drop = FALSE])
	   		rss <- sum(z$residuals^2)
	   		rdf <- z$df.residual
	   	   	se <- sqrt(diag(R) * rss / rdf)
	   	   	est <- z$coefficients
			tval <- est / se
			pval <- 2*pt(abs(tval), rdf, lower.tail = FALSE)
			data.frame(
				q = est[1],
				H = est[2],
				pval = pval[2],
				rsquared = 1 - rss / sum((y - mean(y))^2)
				)
		} else {
			data.frame(
				q = NA,
				H = NA,
				pval = NA,
				rsquared = NA
				)
		}
	}

	.dolm <- function(df){
		lmfit <- lm(logsd ~ logscale, data = df)
		data.frame(
			q = lmfit$coefficients[1],
			H = lmfit$coefficients[2],
			pval = summary(lmfit)$coefficients[2,4],
			rsquared = summary(lmfit)$r.squared
		)
	}

	.setValues <- function(r, values){
		terra::values(r) <- values
		return(r)
	}

	melted_df <- raster_values %>% dplyr::group_by(id) %>% dplyr::group_modify(~ .dofastlm(.)) %>% reshape2::melt(id.vars = "id")

	Hurst_rasters <- melted_df %>% 
		dplyr::group_by(variable) %>% 
		dplyr::group_map(~ .setValues(L_raster, .x$value))
	Hurst_rasters <- do.call(c, Hurst_rasters)
	return(Hurst_rasters)
}