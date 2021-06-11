#' Wrapper function to derive the anisotropy exponent from a non-rotated raster
#' @param rstr `Raster` object
#' @param raster_resolution `numeric` the resolution of `rstr` in meters
#' @param angle_step `numeric`, angular step
#' @param niter `numeric`, number of random rotations
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return a `data.frame`
#' @export
#' @keywords zeta
get_radial_angle <- function(rstr, raster_resolution, angle_step, niter = 64){
	angles <- seq(0, 90 - angle_step, angle_step)
	random_angles <- sample(seq(360), niter)
	random_res <- foreach(random_angle = random_angles, .combine = rbind) %do% {
		.rstr <- rstr %>% rotate_raster(-random_angle) %>% raster::raster() 
		res <- foreach(rotation_angle = angles, .combine = cbind) %do% {
			rotated_raster <- rotate_raster(.rstr, -rotation_angle) # rotating counter-clockwise of theta > 0
			mid_row <- nrow(rotated_raster) %/% 2
			mid_col <- ncol(rotated_raster) %/% 2
			hhcf_x <- get_hhcf_(rotated_raster[c(mid_row - 1, mid_row, mid_row + 1),], raster_resolution, margin = 1, average = TRUE)
			hhcf_y <- get_hhcf_(rotated_raster[, c(mid_col - 1, mid_col, mid_col + 1)], raster_resolution, margin = 2, average = TRUE)
			theta_x <- rotation_angle
			theta_y <- (theta_x + 90)
			alpha_x <- get_all_alpha_(hhcf_x, raster_resolution) 
			alpha_y <- get_all_alpha_(hhcf_y, raster_resolution) 
			colnames(alpha_x) <- paste0(colnames(alpha_x), ".x")
			colnames(alpha_y) <- paste0(colnames(alpha_y), ".y")
			res <- dplyr::bind_cols(alpha_x, alpha_y)
			res <- res %>%
				dplyr::select(
					.data$alpha1.x,
					.data$alpha1.y
				)
			colnames(res) <- gsub(".x", paste0(".", theta_x), colnames(res))
			colnames(res) <- gsub(".y", paste0(".", theta_y), colnames(res))
			return(res)
		}
		alpha1 <- res %>% dplyr::select(dplyr::contains("alpha1.")) %>% as.matrix() %>% c()
		if(all(is.na(alpha1))) {
			return(
				data.frame(alpha1_median = NA, alpha1_mad = NA, alpha1_mean = NA, alpha1_sd = NA, theta_perp = NA)
			)
		}
		alpha1_median <- stats::median(alpha1, na.rm = TRUE)
		alpha1_mean <- mean(alpha1, na.rm = TRUE)
		alpha1_mad <- stats::mad(alpha1, constant = 1, na.rm = TRUE)
		alpha1_sd <- stats::sd(alpha1, na.rm = TRUE)
		theta <- res %>% dplyr::select(dplyr::contains("alpha1.")) %>% colnames() %>% gsub("alpha1.", "", .) %>% as.numeric()
		theta <- theta[!is.na(alpha1)]
		alpha1 <- alpha1[!is.na(alpha1)]
		if(length(alpha1) < 2){
			return(
				data.frame(alpha1_median = NA, alpha1_mad = NA, alpha1_mean = NA, alpha1_sd = NA, theta_perp = NA)
			)
		} 
		ww <- scales::rescale(alpha1)
		ww <- ww / sum(ww)
		dens <- density(theta, weights = ww, bw = "SJ-dpi", kernel = "gaussian", n = 512)
		theta_perp <- (- dens$x[which.max(dens$y)]) %% 180
		(theta_perp <- (theta_perp - random_angle) %% 180)
		data.frame(
			alpha1_median = alpha1_median,
			alpha1_mad = alpha1_mad,
			alpha1_mean = alpha1_mean,
			alpha1_sd = alpha1_sd,
			theta_perp = theta_perp
		)
	}
	res <- random_res %>% 
		na.omit() %>% 
		dplyr::summarise(dplyr::across(dplyr::everything(), modeest::parzen, bw = "SJ-dpi", kernel = "gaussian")) %>% 
		as.list()
	return(res)
}
