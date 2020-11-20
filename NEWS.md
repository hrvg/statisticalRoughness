# statisticalRoughness 0.2

New features

- The following functions were added to allow a 2D Fourier analysis of topography:
	+ `detrend_dem()`
	+ `pad_mat()`
	+ `Hann2D()`
	+ `fft2D()`
- The following functions were added to rotate raster according to a 2D Fourier analysis:
	+ `bin()`
	+ `get_normalized_spectral_power_matrix()`
	+ `filter_spectral_power_matrix()`
	+ `get_fourier_angle()`
	+ `rotate_raster()`

Enhancements

- re-organized the Reference and the Articles
- added one vignette describing the basic use of Fourier analysis functions
- added one vignette describing how to rotate rasters along 2DFT main components

# statisticalRoughness 0.1

Initial release.

- Contains the following functions
	+ `all_rasters_to_polygons()`
	+ `compute_Hurst_rasters()`
	+ `compute_Hurst_rasters_internal()`
	+ `get_mean_res()`
	+ `get_mean_sd_rasters()`
	+ `merge_Hurst_rasters()`
	+ `get_all_R_L()`
	+ `get_factors()`
	